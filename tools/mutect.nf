#!/usr/bin/env nextflow

// STEP 1 - pileups. Do this for the normal and abnormal bam
// Create channel for Pileups? have to run script for tumor and normal
process GETPILEUPSUMMARIES {
    maxForks 3
    container "${params.container_gatk}"

    publishDir "${params.outdir}/intermediates", mode: 'copy'
    
    input:
    path bam_file
    path bai_file
    path small_exac_common_3.hg38.vcf.gz

    output:
    path ("*${params.tumor_namepattern}*.pileup.table"), emit: tumor, optional: true
    path ("*${params.normal_namepattern}*.pileup.table"), emit: normal, optional: true

    script:
    """
    gatk GetPileupSummaries \
    --input ${bam_file} \
    --output ${bam_file.baseName}.pileup.table \
    --variant ${params.ref}/small_exac_common_3.hg38.vcf.gz \
    --intervals ${params.ref}/small_exac_common_3.hg38.vcf.gz \
    --reference ${params.ref}/core_ref_GRCh38_hla_decoy_ebv_modified/genome.fa
    """
}

// STEP 2 - contamination - this uses the outputs from the normal and abnormal pileup
process CALCULATECONTAMINATION {
    container "${params.container_gatk}"

    publishDir "${params.outdir}/vcfs", mode: 'copy'

    input:
    path tumor_pileups_table
    path normal_pileups_table

    output:
    path ("${tumor_pileups_table.simpleName}_contamination_table"), emit: contamination
    path ("${tumor_pileups_table.simpleName}_segmentation_table"), emit: segment

    script:
    """
    gatk CalculateContamination \
    --matched-normal ${normal_pileups_table} \
    --input ${tumor_pileups_table} \
    --output ${tumor_pileups_table.simpleName}_contamination_table \
    --tumor-segmentation ${tumor_pileups_table.simpleName}_segmentation_table
    """
}

// STEP 3 - run MUTECT2. Do without PON
// This runs mutect on a region. The format of the 'region' variable is "chr1:1234567-2345678". Uses interval list of 792 regions
// **** Make a channel for regions. Limit the amount of jobs submitted. *****
process MUTECT2 {

    maxForks 12 //set this when running on local scratch to parallelize; up to the max number of cpus available minus 1
    cpus 1 // set cpu to 1; gatk discourages multithreading
    container "${params.container_gatk}"

    input:
    path tumor_bam
    path tumor_bam_sorted_bai // only necessary when nextflow can't resolve path from symlink
    path normal_bam
    path normal_bam_sorted_bai // only necessary when nextflow can't resolve path from symlink
    each chrom // repeat this process for each item in the chrom channel
    val sample_id 
    val normal_command
    path mutect_idx // nextflow can't resolve the rest of these files from symlink to mutect_idx, input paths to each
    path mutect_idx_fai
    path mutect_idx_dict
    
    output:
    path "${t_n}.mutect", emit: mutect
    path "${output}_${idx}.vcf", emit: vcf

    path "${sample_id}_${chrom}_unfiltered.vcf", emit: vcf
    path "${sample_id}_${chrom}_f1r2.tar.gz", emit: f1r2
    path "${sample_id}_${chrom}_unfiltered.vcf.stats", emit: stats
    path "${sample_id}_${chrom}_unfiltered.vcf.idx", emit: index

    script:
    // normal name (sample_id) is whatever the SM: label in your @RG is in the bam file
    """
    gatk Mutect2 \
    --reference ${params.ref}/core_ref_GRCh38_hla_decoy_ebv_modified/genome.fa \
    -I ${tumor_bam.join(' -I ')} \
    -I ${normal_bam.join(' -I ')} \
    -normal ${normalx} \
    -tumor ${tumourx} \
    -L ${chrom} \
    --annotation FisherStrand \
    --output ${sample_id}_${chrom}_unfiltered.vcf \
    --germline-resource ${params.ref}/af-only-gnomad.hg38.vcf.gz \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter  \
    --f1r2-tar-gz ${sample_id}_${chrom}_f1r2.tar.gz \
    --f1r2-max-depth 600 \
    --dont-use-soft-clipped-bases true
    """
}

// Step 4 - merge vcfs
// Merge the output mutect vcf files into a single file - need to make a list of inputs
// get list of input vcfs
// ********* SAME AS BGZIP/PREPARE VCF??????? ***********
process MERGE {
    container "${params.container_gatk}"

    input: 

    output:

    script:
    """
    idx=0
    # make sure file is empty
    >"input_vcf.txt"
    while read region; do
    input_vcf="${output}_${idx}.vcf"
    echo ${input_vcf} >> input_vcfs.txt
    idx=$((idx + 1))
    done < "${params.ref}/hg38.even.intervals"

    # Initialize the input parameter for MergeVcfs
    inputsVCF=""

    # Loop through each line in the file and add the -I flag for each VCF
    while IFS= read -r line
    do
        inputVCF+="-I $line "

    done < "input_vcfs.txt"

    # Output merged VCF
    mergedvcf="${t_n}.mutect.vcf"

    # Run GATK MergeVcfs with all input VCFs
    gatk MergeVcfs $inputsVCF --OUTPUT ${mergedvcf} --SEQUENCE_DICTIONARY "${dict}"
    """
}

// STEP 5 - merge stats
// Merge the stats output by each individual run of mutect -  use the intervals file as a counter
// get list of input stats. If not in workdir this becomes a too many character input and does not work
process MERGESTATS {
    container "${params.container_gatk}"

    publishDir "${params.outdir}/intermediates", mode: 'copy'

    input: 
    path stats
    val sample_id

    output: 
    path "${sample_id}_unfiltered.vcf.all.stats"

    script:
    """
    echo ${stats.join(' ')}
    gatk MergeMutectStats -stats ${stats.join(' -stats ')} -O ${sample_id}_unfiltered.vcf.all.stats
    """
}

// step 6 - read orientation - just use the intervals file as a counter
// get list of read orientations  -  use the intervals file as a counter
process LEARNORIENTATION {
    container "${params.container_gatk}"

    publishDir "${params.outdir}/intermediates", mode: 'copy'

    input: 
    path f1r2file
    val sample_id

    output: 
    path("${sample_id}_read-orientation-model.tar.gz")

    script:
    """
    gatk LearnReadOrientationModel -I ${f1r2file.join(' -I ')} -O ${sample_id}_read-orientation-model.tar.gz 
    """
    
}

// STEP 7 - add filters. This STEP adds information to the filter column of the vcf
process FILTERMUTECT {
    container "${params.container_gatk}"

    publishDir "${params.outdir}/intermediates", mode: 'copy'

    input: 
    path unfiltered_vcf
    path unfiltered_vcf_index
    path mutect_idx
    path mutect_idx_fai
    path mutect_dict
    path vcf_stats
    path read_orientation_model
    path segmentation_table
    path contamination_table
    val sample_id  

    output:
    path("${sample_id}_filtered.vcf")

    script: 
    """
    gatk FilterMutectCalls \
    --variant ${unfiltered_vcf}  \
    --reference ${genome} \
    --output ${sample_id}_mutect.filters.vcf \
    --contamination-table ${contamination_table.join(' --contamination-table ')} \
    --tumor-segmentation ${segmentation_table.join(' --tumor-segmentation ')}  \
    --orientation-bias-artifact-priors ${read_orientation_model} \
    --stats ${vcf_stats}  \
    --filtering-stats ${sample_id}.mutect.filters.txt \
    --min-slippage-length 8  \
    --pcr-slippage-rate 0.1 \
    --max-events-in-region 2
    """
}

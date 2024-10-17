#!/usr/bin/env nextflow

// STEP 1 - pileups
process GetPileupSummaries {
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.sample_id}/mutect/intermediates", mode:'copy'
    
    input:
    path bam
    path bai
    path basename
    // Reference Files
    path small_exac
    path genome_modified

    output:
    path "*.pileup.table", emit: pileup 

    script:
    """
    gatk GetPileupSummaries \
    --input ${bam} \
    --output ${basename}.pileup.table \
    --variant ${small_exac} \
    --intervals ${small_exac} \
    --reference ${genome_modified}
    """
}

// STEP 2 - contamination 
process CalculateContamination {
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.sample_id}/mutect/vcfs", mode:'copy'

    input:
    path normal_pileups_table
    path tumor_pileups_table

    output:
    path "*.contamination_table", emit: contamination
    path "*.segmentation_table", emit: segmentation

    script:
    """
    gatk CalculateContamination \
    --matched-normal ${normal_pileups_table} \
    --input ${tumor_pileups_table} \
    --output ${t_n}.contamination.table \
    --tumor-segmentation ${t_n}.segmentation.table
    """
}

// STEP 3 - run MUTECT2. Do without PON
// This runs mutect on a region. The format of the 'region' variable is "chr1:1234567-2345678". Uses interval list of 792 regions
process MUTECT2 {
    maxForks 12 //set this when running on local scratch to parallelize; up to the max number of cpus available minus 1
    cpus 1 // set cpu to 1; gatk discourages multithreading
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.sample_id}/mutect/results", mode:'copy'

    input:
    each region // repeat this process for each item in the region channel
    each idx    // repeat this process for each item in the index channel
    path normal_bam
    path tumor_bam
    // Reference File
    path genome_modified
    path gnomad
    
    output:
    path "${t_n}.mutect_${idx}.vcf", emit: vcf
    path "${t_n}.mutect_${idx}.vcf.f1r2.tar.gz", emit: f1r2
    path "${t_n}.mutect_${idx}.vcf.stats", emit: stats
    path "${t_n}.mutect_${idx}.vcf.idx", emit: index

    script:
    // normal name (sample_id) is whatever the SM: label in your @RG is in the bam file
    // @RG SM:SA###### (same sample id in filename)

    // options to get normalx and tumourx:
    // 1. use samtools (samtools view -H <bam> | grep "^@RG" | cut -f3 | cut -d ":" -f2 )
    // 2. include sample id in params file (get from bam filename when generating json files)
    """
    t_n="${tumourx}_${normalx}"

    gatk Mutect2 \
    --reference ${genome_modified} \
    -I ${normal_bam} \
    -I ${tumor_bam} \
    -normal ${normalx} \
    -tumor ${tumourx} \
    -L ${region} \
    --annotation FisherStrand \
    --output ${t_n}.mutect_${idx}.vcf \
    --germline-resource ${gnomad} \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter  \
    --f1r2-tar-gz "${t_n}.mutect_${idx}.vcf.f1r2.tar.gz" \
    --f1r2-max-depth 600 \
    --dont-use-soft-clipped-bases true
    """
}

// Step 4 - merge vcfs
// Merge the output mutect vcf files into a single file 
process MERGE {
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.sample_id}/mutect/merged", mode:'copy'


    input: 
    path vcf_channel // all output vcfs from MUTECT2 process
    // Reference File
    path genome_modified

    output:
    path "${t_n}.mutect.vcf", emit: merged_vcf

    script:
    """
    gatk MergeVcfs \
    -I ${vcf_channel.join(' -I ')} \
    --OUTPUT ${t_n}.mutect.vcf \
    --SEQUENCE_DICTIONARY ${genome_modified}/genome.dict
    """
}

// STEP 5 - merge stats
// Merge the stats output by each individual run of mutect 
process MergeStats {
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.sample_id}/mutect/merged", mode:'copy'

    input: 
    path stats

    output: 
    path "${t_n}.mutect.vcf.stats", emit: merged_stats

    script:
    """
    gatk MergeMutectStats -stats ${stats.join(' -stats ')} --output ${t_n}.mutect.vcf.stats
    """
}

// STEP 6 - read orientation
process LearnOrientation {
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.sample_id}/intermediates", mode: 'copy'

    input: 
    path f1r2file

    output: 
    path "${t_n}.mutect.vcf.f1r2.tar.gz", emit: read_orientation

    script:
    """
    gatk LearnReadOrientationModel -I ${f1r2file.join(' -I ')} --output ${t_n}.mutect.vcf.f1r2.tar.gz
    """
    
}

// STEP 7 - add filters. This STEP adds information to the filter column of the vcf
process FilterMutect {
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.sample_id}/intermediates", mode: 'copy'

    input: 
    path merged_vcf
    path genome
    path contam_table
    path segment_table
    path read_orientation
    path merged_stats

    output:
    path "${t_n}.mutect.filters.vcf"
    path "${t_n}.mutect.filters.txt"

    script: 
    """
    gatk FilterMutectCalls \
    --variant ${merged_vcf}  \
    --reference ${genome}/genome.fa \
    --output ${t_n}.mutect.filters.vcf \
    --contamination-table ${contam_table} \
    --tumor-segmentation ${segment_table}  \
    --orientation-bias-artifact-priors ${read_orientation} \
    --stats ${merged_stats} \
    --filtering-stats ${t_n}.mutect.filters.txt \
    --min-slippage-length 8  \
    --pcr-slippage-rate 0.1 \
    --max-events-in-region 2
    """
}
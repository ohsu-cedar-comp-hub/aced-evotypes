#!/usr/bin/env nextflow

// STEP 1 - pileups
process GetPileupSummaries {
    errorStrategy 'ignore'
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.case_id}/mutect", mode:'copy'
    
    input:
    path files
    val normal_file
    val tumor_file
    // Reference Files
    path small_exac
    path small_exac_tbi
    path genome_modified

    output:
    path "*.normal.pileup.table", emit: normal_pileup 
    path "*.tumor.pileup.table", emit: tumor_pileup

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    normal_basename=\$(basename ${normal_file} .bam)
    tumor_basename=\$(basename ${tumor_file} .bam)

    workdir=\$(pwd)
    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    echo "listing items in files: ${files}"
    ls ${files}

    # Normal BAM
    gatk GetPileupSummaries \
    --input ${files}/${normal_file} \
    --output \${normal_basename}.normal.pileup.table \
    --variant ${small_exac} \
    --intervals ${small_exac} \
    --reference ${genome_modified}/genome.fa

    # Tumor BAM
    gatk GetPileupSummaries \
    --input ${files}/${tumor_file} \
    --output \${tumor_basename}.tumor.pileup.table \
    --variant ${small_exac} \
    --intervals ${small_exac} \
    --reference ${genome_modified}/genome.fa
    """
}

// STEP 2 - contamination 
process CalculateContamination {
    errorStrategy 'ignore'
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.case_id}/mutect", mode:'copy'

    input:
    path normal_pileup
    path tumor_pileup
    val normal_SM
    val tumor_SM

    output:
    path "${tumor_SM}_${normal_SM}.contamination.table", emit: contamination
    path "${tumor_SM}_${normal_SM}.segmentation.table", emit: segmentation

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    gatk CalculateContamination \
    --matched-normal ${normal_pileup} \
    --input ${tumor_pileup} \
    --output ${tumor_SM}_${normal_SM}.contamination.table \
    --tumor-segmentation ${tumor_SM}_${normal_SM}.segmentation.table
    """
}

// STEP 3 - run MUTECT2. Do without PON
// This runs mutect on a region. The format of the 'region' variable is "chr1:1234567-2345678". Uses interval list of 793 regions
process MUTECT2 {
    errorStrategy 'ignore'
    maxForks 20 //set this when running on local scratch to parallelize; up to the max number of cpus available minus 1
    cpus 1 // set cpu to 1; gatk discourages multithreading
    container "${params.container_gatk}"

    input:
    val region
    val idx
    path files
    val normal_file
    val normal_SM
    val tumor_file
    val tumor_SM
    // Reference Files
    path genome_modified
    path gnomad
    path gnomad_tbi
    
    output:
    path "${tumor_SM}_${normal_SM}.mutect_${idx}.vcf", emit: vcf
    path "${tumor_SM}_${normal_SM}.mutect_${idx}.vcf.f1r2.tar.gz", emit: f1r2
    path "${tumor_SM}_${normal_SM}.mutect_${idx}.vcf.stats", emit: stats
    path "${tumor_SM}_${normal_SM}.mutect_${idx}.vcf.idx", emit: index

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    workdir=\$(pwd)
    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    echo "listing items in files: ${files}"
    ls ${files}

    gatk Mutect2 \
    --reference ${genome_modified}/genome.fa \
    -I ${files}/${normal_file} \
    -I ${files}/${tumor_file} \
    -normal ${normal_SM} \
    -tumor ${tumor_SM} \
    -L ${region} \
    --annotation FisherStrand \
    --output ${tumor_SM}_${normal_SM}.mutect_${idx}.vcf \
    --germline-resource ${gnomad} \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter  \
    --f1r2-tar-gz "${tumor_SM}_${normal_SM}.mutect_${idx}.vcf.f1r2.tar.gz" \
    --f1r2-max-depth 600 \
    --dont-use-soft-clipped-bases true
    """
}

// Step 4 - merge vcfs
// Merge the output mutect vcf files into a single file 
process MERGE {
    errorStrategy 'ignore'
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.case_id}/mutect/merged", mode:'copy'

    input: 
    path vcf_channel // all output vcfs from MUTECT2 process
    // Reference File
    path genome_modified
    val normal_SM
    val tumor_SM

    output:
    path "${tumor_SM}_${normal_SM}.mutect.vcf", emit: merged_vcf

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    gatk MergeVcfs \
    -I ${vcf_channel.join(' -I ')} \
    --OUTPUT ${tumor_SM}_${normal_SM}.mutect.vcf \
    --SEQUENCE_DICTIONARY ${genome_modified}/genome.dict
    """
}

// STEP 5 - merge stats
// Merge the stats output by each individual run of mutect 
process MergeStats {
    errorStrategy 'ignore'
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.case_id}/mutect/merged", mode:'copy'

    input: 
    path stats
    val normal_SM
    val tumor_SM

    output: 
    path "${tumor_SM}_${normal_SM}.mutect.vcf.stats", emit: merged_stats

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    gatk MergeMutectStats -stats ${stats.join(' -stats ')} --output ${tumor_SM}_${normal_SM}.mutect.vcf.stats
    """
}

// STEP 6 - read orientation
process LearnOrientation {
    errorStrategy 'ignore'
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.case_id}/mutect/merged", mode: 'copy'

    input: 
    path f1r2file
    val normal_SM
    val tumor_SM

    output: 
    path "${tumor_SM}_${normal_SM}.mutect.vcf.f1r2.tar.gz", emit: read_orientation

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    gatk LearnReadOrientationModel -I ${f1r2file.join(' -I ')} --output ${tumor_SM}_${normal_SM}.mutect.vcf.f1r2.tar.gz
    """
    
}

// STEP 7 - add filters. This STEP adds information to the filter column of the vcf
process FilterMutect {
    errorStrategy 'ignore'
    container "${params.container_gatk}"
    publishDir "${params.bucket}/${params.case_id}/mutect/filtered", mode: 'copy'

    input: 
    path merged_vcf
    path genome
    path contam_table
    path segment_table
    path read_orientation
    path merged_stats
    val normal_SM
    val tumor_SM

    output:
    path "${tumor_SM}_${normal_SM}.mutect.filters.vcf"
    path "${tumor_SM}_${normal_SM}.mutect.filters.txt"

    script: 
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    gatk FilterMutectCalls \
    --variant ${merged_vcf}  \
    --reference ${genome}/genome.fa \
    --output ${tumor_SM}_${normal_SM}.mutect.filters.vcf \
    --contamination-table ${contam_table} \
    --tumor-segmentation ${segment_table}  \
    --orientation-bias-artifact-priors ${read_orientation} \
    --stats ${merged_stats} \
    --filtering-stats ${tumor_SM}_${normal_SM}.mutect.filters.txt \
    --min-slippage-length 8  \
    --pcr-slippage-rate 0.1 \
    --max-events-in-region 2
    """
}
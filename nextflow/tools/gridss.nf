#!/usr/bin/env nextflow

// Step 1
process GRIDSS {
    container "${container_gridss}"
    publishDir "${params.bucket}/${params.sample_id}/gridss", mode:'copy'


    input:
    path normal_bam
    path tumor_bam
    // Reference Files
    path genome

    output:
    path "${seq_id_t}_v_${seq_id_n}.gridss.vcf", emit: vcf
    path "${t_n}.gridss.assembly"

    script:
    """
    java -jar $GRIDSS_JAR gridss \
    -r ${genome} \
    -o ${seq_id_t}_v_${seq_id_n}.gridss.vcf \
    -a ${t_n}.gridss.assembly \
    ${normal_bam} ${tumor_bam}
    """
}

process GridssSomaticFilter {
    container "${container_gridss}"
    publishDir "${params.bucket}/${params.sample_id}/gridss", mode:'copy'


    input:
    path vcf
    // Reference file
    path BSgenome

    output:
    path "${seq_id_t}_v_${seq_id_n}.gridss.filtered.vcf"
    path "${seq_id_t}_v_${seq_id_n}.gridss.semifiltered.vcf"

    script:
    """
    {rscript} /script/dir/gridss_somatic_filter \
    --pondir gridss-2.13.2 \
    --ref ${BSgenome}/BSgenome.Hsapiens.UCSC.hg38 \
    --input ${vcf} \
    --output ${seq_id_t}_v_${seq_id_n}.gridss.filtered.vcf \
    --fulloutput ${seq_id_t}_v_${seq_id_n}.gridss.semifiltered.vcf \
    --normalordinal 1 \
    --tumourordinal 2 \
    --s "/script/dir/gridss_trial_17_5" \
    --gc
    """
}
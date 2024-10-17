#!/usr/bin/env nextflow

// Step 1
process GRIDSS {
    container "${params.container_gridss}"
    publishDir "${params.bucket}/${params.sample_id}/gridss", mode:'copy'

    input:
    path normal_bam
    path tumor_bam
    // Reference File
    path genome

    output:
    path "${t_n}.gridss.vcf", emit: vcf
    path "${t_n}.gridss.assembly"

    script:
    """
    workdir=\${pwd}

    gridss \
    --jvmheap "15g" \
    --workingdir \${workdir} \
    --reference ${genome} \
    --assembly \${workdir}/${t_n}.gridss.assembly \
    --output \${workdir}/${t_n}.gridss.vcf \
    ${normal_bam} ${tumor_bam}
    """
}

process GridssSomaticFilter {
    container "${params.container_gridss_r}"
    publishDir "${params.bucket}/${params.sample_id}/gridss/filtered", mode:'copy'

    input:
    path vcf
    path gridss_dir
    path pondir

    output:
    path "${t_n}.gridss.filtered.vcf"
    path "${t_n}.gridss.semifiltered.vcf"

    script:
    """
    {rscript} ${gridss_dir}/gridss_somatic_filter \
    --pondir ${pondir} \
    --ref BSgenome.Hsapiens.UCSC.hg38 \
    --input ${vcf} \
    --output ${t_n}.gridss.filtered.vcf \
    --fulloutput ${t_n}.gridss.semifiltered.vcf \
    --normalordinal 1 \
    --tumourordinal 2 \
    --s ${gridss_dir} \
    --gc
    """
}
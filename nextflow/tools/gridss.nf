#!/usr/bin/env nextflow

// Step 1
process GRIDSS {
    errorStrategy 'ignore'
    container "${params.container_gridss}"
    publishDir "${params.bucket}/${params.case_id}/gridss", mode:'copy'

    input:
    path files
    val normal_file
    val normal_SM
    val tumor_file
    val tumor_SM
    // Reference File
    path core_ref

    output:
    path "${tumor_SM}_${normal_SM}.gridss.vcf", emit: vcf
    path "${tumor_SM}_${normal_SM}.gridss.assembly"

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
    
    gridss \
    --jvmheap "15g" \
    --workingdir \${workdir} \
    --reference ${core_ref}/genome.fa \
    --assembly \${workdir}/${tumor_SM}_${normal_SM}.gridss.assembly \
    --output \${workdir}/${tumor_SM}_${normal_SM}.gridss.vcf \
    -t 8 \
    ${files}/${normal_file} ${files}/${tumor_file}
    """
}

process GridssSomaticFilter {
    errorStrategy 'ignore'
    container "${params.container_gridss_r}"
    publishDir "${params.bucket}/${params.case_id}/gridss/filtered", mode:'copy'

    input:
    path vcf
    path gridss_dir
    path pondir
    val normal_SM
    val tumor_SM

    output:
    path "${tumor_SM}_${normal_SM}.gridss.filtered.vcf"
    path "${tumor_SM}_${normal_SM}.gridss.semifiltered.vcf"

    script:
    """
    Rscript ${gridss_dir}/gridss_somatic_filter \
    --input ${vcf} \
    --pondir ${pondir} \
    --ref BSgenome.Hsapiens.UCSC.hg38 \
    --output ${tumor_SM}_${normal_SM}.gridss.filtered.vcf \
    --fulloutput ${tumor_SM}_${normal_SM}.gridss.semifiltered.vcf \
    --normalordinal 1 \
    --tumourordinal 2 \
    -s ${gridss_dir} \
    --gc
    """
}
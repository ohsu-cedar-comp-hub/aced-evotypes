#!/usr/bin/env nextflow

// run Manta
// Step 1
process ConfigManta {
    errorStrategy 'ignore'
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.case_id}", mode:'copy'

    input: 
    path files
    val normal_file
    val tumor_file
    // Reference File
    path core_ref

    output:
    path "manta/*"
    path "manta", emit: manta

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
    mkdir -p \${workdir}/manta

    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    echo "listing items in files: ${files}"
    ls ${files}

    configManta.py \
    --normalBam ${files}/${normal_file} \
    --tumorBam ${files}/${tumor_file} \
    --referenceFasta ${core_ref}/genome.fa \
    --runDir \${workdir}/manta
    """
}
    
// Step 2
process MantaRunWorkflow {
    errorStrategy 'ignore'
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.case_id}", mode:'copy'

    input:
    path files
    path manta
    // Reference File
    path core_ref

    output:
    path "manta/*.txt"
    path "manta/results/*"
    path "${manta}/results/variants/somaticSV.vcf.gz", emit: somaticSV

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    ${manta}/runWorkflow.py -j 8
    """
}

// SV quality filtering
process MantaQualityFilter {
    errorStrategy 'ignore'
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.case_id}/manta/filtered", mode:'copy'

    input: 
    path somaticSV
    val normal_SM
    val tumor_SM

    output:
    path "MANTAPASS_${tumor_SM}_${normal_SM}.vcf"

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

    bcftools filter \
    -O v  \
    -o \$workdir/MANTAPASS_${tumor_SM}_${normal_SM}.vcf \
    -i "FILTER == 'PASS'" \
    "${somaticSV}"
    """
}
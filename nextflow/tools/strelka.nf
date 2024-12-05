#!/usr/bin/env nextflow

// run Strelka
// Step 1
process ConfigStrelka {
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
    path "strelka/*"
    path "strelka", emit: strelka

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
    mkdir -p \${workdir}/strelka

    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    echo "listing items in files: ${files}"
    ls ${files}

    configureStrelkaSomaticWorkflow.py \
    --normalBam ${files}/${normal_file} \
    --tumorBam ${files}/${tumor_file} \
    --outputCallableRegions \
    --referenceFasta ${core_ref}/genome.fa \
    --runDir \${workdir}/strelka
    """
}

// Step 2 
process StrelkaRunWorkflow {
    errorStrategy 'ignore'
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.case_id}", mode:'copy'

    input: 
    path files
    path strelka
    // Reference file
    path core_ref

    output:
    path "strelka/*.txt"
    path "strelka/results/*"
    path "${strelka}/results/variants/somatic.indels.vcf.gz", emit: somatic_indels
    path "${strelka}/results/variants/somatic.snvs.vcf.gz", emit: somatic_snvs

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    ${strelka}/runWorkflow.py -m local -j 8
    """
}

// INDEL quality filtering
process StrelkaIndelFilter {
    errorStrategy 'ignore'
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.case_id}/strelka/filtered", mode:'copy'

    input:
    path somatic_indels
    val normal_SM
    val tumor_SM

    output:
    path "STRELKAPASS_somatic.indels_${tumor_SM}_${normal_SM}.vcf"

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
    -o \${workdir}/STRELKAPASS_somatic.indels_${tumor_SM}_${normal_SM}.vcf \
    -i "FILTER == 'PASS'" \
    "${somatic_indels}"
    """
}

// SNV quality filtering
process StrelkaSnvFilter {
    errorStrategy 'ignore'
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.case_id}/strelka/filtered", mode:'copy'

    input:
    path somatic_snvs
    val normal_SM
    val tumor_SM

    output:
    path "STRELKAPASS_somatic.snvs_${tumor_SM}_${normal_SM}.vcf"

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
    -o \${workdir}/STRELKAPASS_somatic.snvs_${tumor_SM}_${normal_SM}.vcf \
    -i "FILTER == 'PASS'" \
    "${somatic_snvs}"
    """
}
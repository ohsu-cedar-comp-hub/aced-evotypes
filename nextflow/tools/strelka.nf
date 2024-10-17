#!/usr/bin/env nextflow

// run Strelka
// Step 1
process ConfigStrelka {
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.sample_id}", mode:'copy'

    input:
    path normal_bam
    path normal_bai
    path tumor_bam
    path tumor_bai
    // Reference File
    path core_ref

    output:
    path "strelka/*"
    path "strelka", emit: strelka

    script:
    """
    workdir=\$(pwd)
    mkdir \${workdir}/strelka

    configureStrelkaSomaticWorkflow.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --outputCallableRegions \
    --referenceFasta ${core_ref}/genome.fa \
    --runDir \${workdir}/strelka
    """
}

// Step 2 
process StrelkaRunWorkflow {
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.sample_id}", mode:'copy'

    input: 
    path strelka

    output:
    path "${strelka}/results/*"
    path "${strelka}/results/variants/somatic.indels.vcf.gz", emit: somatic_indels
    path "${strelka}/results/variants/somatic.snvs.vcf.gz", emit: somatic_snvs

    script:
    """
    ${strelka}/runWorkflow.py -m local -j $CPUS
    """
}

// INDEL quality filtering
process StrelkaIndelFilter {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/strelka/filtered", mode:'copy'

    input:
    path somatic_indels

    output:
    path "STRELKAPASS_somatic.indels_*"

    script:
    """
    workdir=\$(pwd)

    bcftools filter \
    -O v  \
    -o \${workdir}/STRELKAPASS_somatic.indels_${t_n}.vcf \
    -i "FILTER == 'PASS'" \
    "${somatic_indels}"
    """
}

// SNV quality filtering
process StrelkaSnvFilter {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/strelka/filtered", mode:'copy'

    input:
    path somatic_snvs

    output:
    path "STRELKAPASS_somatic.snvs_*"

    script:
    """
    workdir=\$(pwd)

    bcftools filter \
    -O v  \
    -o \${workdir}/STRELKAPASS_somatic.snvs_${t_n}.vcf \
    -i "FILTER == 'PASS'" \
    "${somatic_snvs}"
    """
}
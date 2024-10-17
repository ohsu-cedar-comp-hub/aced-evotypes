#!/usr/bin/env nextflow

// run Manta
// Step 1
process ConfigManta {
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
    path "manta/*"
    path "manta", emit: manta

    script:
    """
    workdir=\$(pwd)
    mkdir \${workdir}/manta

    configManta.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --referenceFasta ${core_ref}/genome.fa \
    --runDir \${workdir}/manta
    """
}

// Step 2
process MantaRunWorkflow {
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.sample_id}", mode:'copy'

    input:
    path manta

    output:
    path "${manta}/results/*"
    path "${manta}/results/variants/somaticSV.vcf.gz", emit: somaticSV


    script:
    """
    ${manta}/runWorkflow.py -j $CPUS
    """
}

// SV quality filtering
process MantaQualityFilter {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/manta/filtered", mode:'copy'

    input: 
    path somaticSV

    output:
    path "MANTAPASS_${t_n}.vcf"

    script:
    """
    workdir=\$(pwd)

    bcftools filter \
    -O v  \
    -o \$workdir/MANTAPASS_${t_n}.vcf \
    -i "FILTER == 'PASS'" \
    "${somaticSV}"
    """
}
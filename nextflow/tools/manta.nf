#!/usr/bin/env nextflow

// run Manta
// Step 1
process ConfigManta {
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.sample_id}/manta", mode:'copy'

    input: 
    path normal_bam
    path tumor_bam
    // Reference File
    path genome

    output:
    path "manta_results", emit: manta_results //this is a directory

    script:
    """
    workdir=\$(pwd)

    configManta.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --referenceFASTA ${genome} \
    --runDir \${workdir}/manta_results
    """
}

// Step 2
process MantaRunWorkflow {
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.sample_id}/manta", mode:'copy'


    input:
    path manta_results

    output:
    path "results", emit: results // results contains 3 directories (evidence, stats and variants)

    script:
    """
    ${manta_results}/runWorkflow.py -j $CPUS
    """
}

// SV quality filtering
process MantaQualityFilter {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/manta", mode:'copy'

    input: 
    path results

    output:
    MANTAPASS_${t_n}.vcf

    script:
    """
    workdir=\$(pwd)

    bcftools filter \
    -O v  \
    -o \${workdir}/MANTAPASS_${t_n}.vcf \
    -i "FILTER == 'PASS'" \
    "${results}/variants/somaticSV.vcf.gz"
    """
}
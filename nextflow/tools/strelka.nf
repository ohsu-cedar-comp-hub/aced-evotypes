#!/usr/bin/env nextflow

// run Strelka
// Step 1
process ConfigStrelka {
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.sample_id}/strelka", mode:'copy'

    input:
    path normal_bam
    path tumor_bam
    // Reference File
    path genome

    output:
    path "strelka_results", emit: strelka_results

    script:
    """
    workdir=\$(pwd)

    configureStrelkaSomaticWorkflow.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --outputCallableRegions \
    --referenceFasta ${genome} \
    --runDir \${workdir}/strelka_results
    """
}

// Step 2 
process StrelkaRunWorkflow {
    container "${params.container_strelka_manta}"
    publishDir "${params.bucket}/${params.sample_id}/strelka", mode:'copy'


    input: 
    path strelka_results

    output:
    path "results", emit: results // results contains 2 directories (stats and variants)

    script:
    """
    /strelka_results/runWorkflow.py -m local -j $CPUS
    """
}

// INDEL quality filtering
process StrelkaIndelFilter {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/strelka", mode:'copy'

    input:
    path results

    output:
    path "STRELKAPASS_somatic.indels_${t_n}.vcf"

    script:
    """
    workdir=\$(pwd)

    bcftools filter \
    -O v  \
    -o "\${workdir}/STRELKAPASS_somatic.indels_${t_n}.vcf" \
    -i "FILTER == 'PASS'" \
    "${results}/variants/somatic.indels.vcf.gz"
    """
}

// SNV quality filtering
process StrelkaSnvFilter {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/strelka", mode:'copy'

    input:
    path results

    output:
    path "STRELKAPASS_somatic.snvs_${t_n}.vcf"

    script:
    """
    workdir=\$(pwd)

    bcftools filter \
    -O v  \
    -o "\${workdir}/STRELKAPASS_somatic.snvs_${t_n}.vcf" \
    -i "FILTER == 'PASS'" \
    "${results}/variants/somatic.snvs.vcf.gz"
    """
}
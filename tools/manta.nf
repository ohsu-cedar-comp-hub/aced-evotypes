#!/usr/bin/env nextflow

// run Manta
// Step 1
process ConfigManta {
    container "${params.container_strelka_manta}"

    input: 
    path normal_cram_realign
    path tumor_cram_realign

    output:
    path /manta_results, emit: manta_results

    script:
    """
    configManta.py \
    --normalBam /data/<normal cram> \
    --tumorBam /data/<tumor cram> \
    --referenceFASTA ${params.ref}/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
    --runDir /manta_results

    """
}

// Step 2
process MantaRunWorkflow {
    container "${params.container_strelka_manta}"

    input:
    path ConfigManta.out.manta_results

    output:
    manta output files

    script:
    """
    /manta_results/runWorkflow.py -j $CPUS
    """
}

// SV quality filtering
process MantaQualityFilter {
    container "${params.container_bcftools}"

    input: 
    /results/somaticSV.vcf.gz

    output:
    /results/MANTAPASS_${t_n}.vcf

    script:
    """
    bcftools filter \
    -O v  \
    -o ${workplace}/MANTAPASS_${t_n}.vcf \
    -i "FILTER == 'PASS'" \
    "${workplace}/results/variants/somaticSV.vcf.gz"
    """
}
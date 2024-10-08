#!/usr/bin/env nextflow

// run Strelka
// Step 1
process ConfigStrelka {
    container "${params.container_strelka_manta}"
    
    input:
    normal bam
    tumor bam
    /ref/genome.fa

    output:
    /strelka_results

    script:
    """
    configureStrelkaSomaticWorkflow.py \
    --normalBam /data/<normal bam> \
    --tumorBam /data/<tumor bam> \
    --outputCallableRegions \  # is this needed?
    --referenceFasta /ref/genome.fa \  # core_re_GRCh38_hla_decoy_ebv/genome.fa
    --runDir /strelka_results
    """
}

// Step 2 
process StrelkaRunWorkflow {
    container "${params.container_strelka_manta}"

    input: 
    /strelka_results

    output:
    /strelka_results

    script:
    """
    /strelka_results/runWorkflow.py -m local -j $CPUS
    """
}

// INDEL quality filtering
process StrelkaIndelFilter {
    container "${params.container_bcftools}"

    input:
    /strelka_results/somatic.indels.vcf.gz

    output:
    STRELKAPASS_somatic.indels_${t_n}.vcf

    script:
    """
    bcftools filter \
    -O v  \
    -o "${workplace}/STRELKAPASS_somatic.indels_${t_n}.vcf" \
    -i "FILTER == 'PASS'" \
    "${workplace}/results/variants/somatic.indels.vcf.gz"
    """
}

// SNV quality filtering
process StrelkaSnvFilter {
    container "${params.container_bcftools}"
    
    input:
    /results/variants/somatic.snvs.vcf.gz

    output:
    STRELKAPASS_somatic.snvs_${t_n}.vcf

    script:
    """
    bcftools filter \
    -O v  \
    -o "${workplace}/STRELKAPASS_somatic.snvs_${t_n}.vcf" \
    -i "FILTER == 'PASS'" \
    "${workplace}/results/variants/somatic.snvs.vcf.gz"
    """
}
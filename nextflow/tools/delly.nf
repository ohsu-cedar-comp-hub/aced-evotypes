#!/usr/bin/env nextflow

// run allele count steps
process DellyCall {
    container "${params.container_delly}"
    publishDir "${params.bucket}/${params.sample_id}/delly", mode:'copy'

    input: 
    path normal_bam
    path normal_bai
    path tumor_bam
    path tumor_bai
    // Reference Files
    path genome
    path hg38_excl

    output:
    path "*.bcf", emit: bcf

    script:
    """
    delly call \
    -g ${genome} \
    -o delly_${t_n}.bcf \
    -x ${hg38_excl} \
    ${tumor_bam} \
    ${normal_bam}
    """
}

// file has both germline and somatic and must be filtered
// first need to make a table with the headers from the vcf file for the tumor and control
process DellyTable {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/delly", mode:'copy'

    input:
    path bcf

    output:
    path "samples.txt", emit: samples 

    script:
    """
    bcftools view ${bcf} | grep '^#CHROM' | awk -F '\t' '{print $10 "\ttumor\n"$11"\tcontrol"}' > samples.txt 
    """
}

// now can run Delly's optimised filtering that gets rid of any with possible evidence in controls
process DellyFilter {
    container "${params.container_delly}"
    publishDir "${params.bucket}/${params.sample_id}/delly", mode:'copy'

    input: 
    path bcf
    path samples

    output: 
    path "Bdelly_${t_n}.bcf", emit: Bdelly

    script:
    """
    delly filter \
    -f somatic \
    -s ${samples} \
    -o Bdelly_${t_n}.bcf \
    ${bcf}
    """
}

// SV quality filtering
process DellyQualityFilter {
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.sample_id}/delly", mode:'copy'

    input:
    path Bdelly

    output:
    path "dellyPASS_${t_n}.vcf"

    script:
    """
    bcftools filter \
    -O v \
    -o dellyPASS_${t_n}.vcf \
    -i "FILTER == 'PASS'" \
    ${Bdelly}
    """
}
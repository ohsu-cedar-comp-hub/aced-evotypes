#!/usr/bin/env nextflow

//Run the allele_count steps
process DellyCall {
    container "${params.container_delly}"

    input: 
    path normal_cram_realign
    path tumor_cram_realign

    output:
    path "delly_${t_n}.bcf", emit: bcf

    script:
    """
    delly call \
    -g ${params.ref}/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
    -o /data/output/delly_${t_n}.bcf \
    -x ${params.ref}/human.hg38.excl.tsv \
    /data/tmdata/${tumorx}.bam \
    /data/nmdata/${normalx}.bam
    """
}
//file has both germline and somatic and must be filtered

//first need to make a table with the headers from the vcf file for the tumor and control
process DellyTable {
    container "${params.container_bcftools}"

    input:
    path DellyCall.out.bcf

    output:
    path "${workplace}/samples.txt", emit: samples 

    script:
    """
    bcftools view ${workplace}/delly_${t_n}.bcf|grep '^#CHROM' | awk -F '\t' '{print $10 "\ttumor\n"$11"\tcontrol"}'>${workplace}/samples.txt 
    """
}

// now can run Delly's optimised filtering that gets rid of any with possible evidence in controls
process DellyFilter {
    container "${params.container_delly}"

    input: 
    path DellyTable.out.samples
    path DellyCall.out.bcf

    output: 
    path "Bdelly_${t_n}.bcf", emit: Bdelly

    script:
    """
    delly filter -f somatic -s /data/output/samples.txt  -o /data/output/Bdelly_${t_n}.bcf /data/output/delly_${t_n}.bcf
    """
}

// SV quality filtering
process DellyQualityFilter {
    container "${params.container_bcftools}"
    
    input:
    path DellyFilter.out.Bdelly

    output:
    path "${workplace}/dellyPASS_${t_n}.vcf", emit: dellyPass 

    script:
    """
    conda activate bcftools
    
    bcftools filter \
    -O v \
    -o ${workplace}/dellyPASS_${t_n}.vcf \
    -i "FILTER == 'PASS'" \
    ${workplace}/Bdelly_${t_n}.bcf
    """
}
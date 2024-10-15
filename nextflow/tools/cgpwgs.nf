#!/usr/bin/env nextflow

process CGPWGS {
    container "${params.container_cgpwgs}"
    publishDir "${params.bucket}/${params.sample_id}/cgpwgs", mode:'copy'

    input: 
    path normal_bam
    path normal_bai
    path tumor_bam
    path tumor_bai
    // Reference Files
    path core_zip
    path vagrent_zip
    path snv_zip
    path cnv_zip
    path qc_zip

    output:
    path "*.result.tar.gz" // includes ascat/* ; brass/* ; caveman/* ; genotyped/* ; pindel/*

    script:
    """
    workdir=\$(pwd)

    ds-cgpwgs.pl \
    -c 32 \
    -r ${core_zip} \
    -a ${vagrent_zip} \
    -si ${snv_zip} \
    -cs ${cnv_zip} \
    -qc ${qc_zip} \
    -e 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV' \
    -t ${tumor_bam} \
    -tidx ${tumor_bai} \
    -n ${normal_bam} \
    -nidx ${normal_bai} \
    -o \${workdir}
    """
}
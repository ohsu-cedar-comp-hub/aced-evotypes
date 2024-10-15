#!/usr/bin/env nextflow

process CGPMAP {
    container "${params.container_cgpmap}"
    publishDir "${params.bucket}/${params.sample_id}/bam", mode:'copy'

    input: 
    val filename
    val basename
    path cram
    path core_zip
    path bwa_zip

    output:
    path "*realign.OHSU.bam*"
    path "*realign.OHSU.bam", emit: bam 
    path "*flagstat.txt"

    script:
    """
    workdir=\$(pwd)
   
    ds-cgpmap.pl \
    -r ${core_zip} \
    -i ${bwa_zip} \
    -s ${basename}.realign.OHSU \
    -o \${workdir} \
    -t 6 \
    -bwakit \
    -bm2 \
    ${cram}

    samtools flagstat ${cram} > ${basename}.flagstat.txt
    samtools flagstat ${basename}.realign.OHSU > ${basename}.realign.OHSU.flagstat.txt
    """
}
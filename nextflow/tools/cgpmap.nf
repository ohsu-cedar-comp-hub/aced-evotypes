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
    path "bam*"
    path "*.bam", emit: bam 
    path "*flagstat.txt"

    script:
    """
    workdir=\$(pwd)
   
    ds-cgpmap.pl \
    -r ${core_zip} \
    -i ${bwa_zip} \
    -s ${basename} \
    -o \${workdir} \
    -t 6 \
    -bwakit \
    -bm2 \
    ${cram}

    echo "Generating stats for input cram: ${cram}"
    samtools flagstat ${cram} > ${basename}.flagstat.txt

    echo "Generating stats for realigned bam: ${basename}.bam"
    samtools flagstat ${basename}.bam > ${basename}.realign.flagstat.txt
    """
}
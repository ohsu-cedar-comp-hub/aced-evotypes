#!/usr/bin/env nextflow

process CGPMAP {
    container "${params.container_cgpmap}"
    publishDir "${params.bucket}/${params.case_id}/bam", mode:'copy'

    input: 
    val basename
    path cram
    path core_zip
    path bwa_zip

    output:
    path "*bam*"
    path "${basename}.bam", emit: bam 
    path "*flagstat.txt"


    script:
    """
    workdir=\$(pwd)
    echo "workdir: \${workdir}"
   
    echo "Listing files in workdir before analysis"
    ls \${workdir}
    echo ""
    
    SM="\$(samtools view -H \$workdir/${cram} | grep '^@RG' | head -1 | grep -oP 'SM:\\K[^\\s]+')"

    ds-cgpmap.pl \
    -r ${core_zip} \
    -i ${bwa_zip} \
    -s \${SM} \
    -o \${workdir} \
    -t 16 \
    -bwakit \
    -bm2 \
    ${cram}

    echo "Listing files in workdir after alignment"
    ls \${workdir}
    echo ""

    mv \$workdir/\${SM}.bam \$workdir/${basename}.bam
    mv \$workdir/\${SM}.bam.bai \$workdir/${basename}.bam.bai
    mv \$workdir/\${SM}.bam.bas \$workdir/${basename}.bam.bas
    mv \$workdir/\${SM}.bam.maptime \$workdir/${basename}.bam.maptime
    mv \$workdir/\${SM}.bam.md5 \$workdir/${basename}.bam.md5
    mv \$workdir/\${SM}.bam.met \$workdir/${basename}.bam.met

    echo "Listing files in workdir after renaming output"
    ls \${workdir}
    echo ""

    echo "Generating stats for input cram: \$workdir/${cram}"
    samtools flagstat \$workdir/${cram} > \$workdir/${basename}.flagstat.txt

    echo "Generating stats for realigned bam: \$workdir/${basename}.bam"
    samtools flagstat \$workdir/${basename}.bam > \$workdir/${basename}.realign.flagstat.txt
    """
}



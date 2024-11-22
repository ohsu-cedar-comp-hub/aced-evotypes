#!/usr/bin/env nextflow

process CGPWGS {
    container "${params.container_cgpwgs}"
    containerOptions "--home \$(pwd)/cgpwgs"
    publishDir "${params.bucket}/${params.case_id}/cgpwgs", mode:'copy'

    input: 
    path files
    val normal_file
    val tumor_file
    // Reference Files
    path core_zip
    path vagrent_zip
    path snv_zip
    path cnv_zip
    path qc_zip


    output:
    path "cgpwgs/WGS*"
    path "cgpwgs/flag*"
    path "cgpwgs/run.params"
    path "cgpwgs/timings*"

    script:
    """
    workdir=\$(pwd)
    mkdir -p \${workdir}/cgpwgs

    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    echo "listing items in ${files}"
    ls ${files}
    echo ""

    cd \${workdir}/cgpwgs

    ds-cgpwgs.pl \
    -c 30 \
    -r \${workdir}/${core_zip} \
    -a \${workdir}/${vagrent_zip} \
    -si \${workdir}/${snv_zip} \
    -cs \${workdir}/${cnv_zip} \
    -qc \${workdir}/${qc_zip} \
    -e 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV' \
    -t \${workdir}/${files}/${tumor_file} \
    -tidx \${workdir}/${files}/${tumor_file}.bai \
    -n \${workdir}/${files}/${normal_file} \
    -nidx \${workdir}/${files}/${normal_file}.bai \
    -o \${workdir}/cgpwgs

    """
}
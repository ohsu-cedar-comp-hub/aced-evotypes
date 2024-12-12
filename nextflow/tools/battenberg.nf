#!/usr/bin/env nextflow

process AlleleCounter {
    errorStrategy 'ignore'
    container "${params.container_allelecount}"
    publishDir "${params.bucket}/${params.case_id}/battenberg", mode:'copy'

    input:
    path files
    path bin
    path g1000alleles
    val normal_SM
    val tumor_SM

    output:
    path "allele_counting/*"

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    workdir=\$(pwd)
    mkdir -p \${workdir}/allele_counting
    
    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    echo "listing items in files: ${files}"
    ls ${files}
    echo ""

    Rscript ${bin}/runAC_forBattenberg.R \
    -t ${tumor_SM} \
    -n ${normal_SM} \
    --nb "\${workdir}/${files}/${normal_SM}.bam" \
    --tb "\${workdir}/${files}/${tumor_SM}.bam" \
    --cpu 4 \
    -o "\${workdir}/allele_counting" \
    -p "\${workdir}/${g1000alleles}/1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"


    echo ""
    echo "listing workdir after script completion: \${workdir}"
    ls \$workdir
    echo ""
    """
}

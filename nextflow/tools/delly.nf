#!/usr/bin/env nextflow

// run allele count steps
process DellyCall {
    errorStrategy 'ignore'
    container "${params.container_delly}"
    publishDir "${params.bucket}/${params.case_id}", mode:'copy'

    input: 
    path files
    val normal_file
    val normal_SM
    val tumor_file
    val tumor_SM
    // Reference Files
    path genome
    path hg38_excl

    output:
    path "delly/*"
    path "delly", emit: delly
    path "delly/delly_${tumor_SM}_${normal_SM}.bcf", emit: bcf

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
    mkdir -p \${workdir}/delly
    
    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    echo "listing items in files: ${files}"
    ls ${files}

    delly call \
    -g ${genome} \
    -o \${workdir}/delly/delly_${tumor_SM}_${normal_SM}.bcf \
    -x ${hg38_excl} \
    ${files}/${tumor_file} \
    ${files}/${normal_file}

    """
}

// file has both germline and somatic and must be filtered
// first need to make a table with the headers from the vcf file for the tumor and control
process DellyTable {
    errorStrategy 'ignore'
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.case_id}/delly", mode:'copy'

    input:
    path bcf

    output:
    path "samples.txt", emit: samples 

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    bcftools view ${bcf} | grep '^#CHROM' | awk -F '\\t' '{print \$10 "\\ttumor\\n" \$11 "\\tcontrol"}' > samples.txt 
    """
} 

// now can run Delly's optimised filtering that gets rid of any with possible evidence in controls
process DellyFilter {
    errorStrategy 'ignore'
    container "${params.container_delly}"
    publishDir "${params.bucket}/${params.case_id}", mode:'copy'

    input: 
    path delly
    path samples
    val normal_SM
    val tumor_SM

    output: 
    path "delly/filtered/*"
    path "delly/filtered/Bdelly_${tumor_SM}_${normal_SM}.bcf", emit: Bdelly

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
    mkdir -p \${workdir}/delly/filtered
    
    echo "listing workdir: \${workdir}"
    ls \$workdir
    echo ""

    delly filter \
    -f somatic \
    -s ${samples} \
    -o \${workdir}/delly/filtered/Bdelly_${tumor_SM}_${normal_SM}.bcf \
    ${delly}/delly_${tumor_SM}_${normal_SM}.bcf
    """
}

// SV quality filtering
process DellyQualityFilter {
    errorStrategy 'ignore'
    container "${params.container_bcftools}"
    publishDir "${params.bucket}/${params.case_id}/delly/filtered", mode:'copy'

    input:
    path Bdelly
    val normal_SM
    val tumor_SM

    output:
    path "dellyPASS_${tumor_SM}_${normal_SM}.vcf"

    script:
    """
    echo "NEXTFLOW PIPELINE VERSION"
    echo "****************************"
    echo "params.release: ${params.release}"
    echo "params.releasedate: ${params.releasedate}"
    echo "params.githublink: ${params.githublink}"
    echo "****************************"
    echo ""

    bcftools filter \
    -O v \
    -o dellyPASS_${tumor_SM}_${normal_SM}.vcf \
    -i "FILTER == 'PASS'" \
    ${Bdelly}
    """
}
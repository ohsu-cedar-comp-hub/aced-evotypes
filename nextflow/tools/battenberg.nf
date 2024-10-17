#!/usr/bin/env nextflow

process AlleleCounter {
    // load R env with required packages
    container "${params.container_ac_r}"
    publishDir "${params.bucket}/${params.sample_id}/battenberg", mode:'copy'

    input:
    path normal_bam
    path tumor_bam

    output:
    path "allele_counting/*"

    script:
    """
    workdir=\${pwd}
    mkdir \${workdir}/allele_counting

    # execute alleleCounting
    R CMD BATCH "--no-restore-data --no-save --args \
    -t ${tumour_id} \
    -n ${normal_id} \
    --nb ${normal_bam} \
    --tb ${tumour_bam} \
    --cpu 16 \
    -o \${workdir}/allele_counting "
    runAC_forBattenberg.R "\${workdir}/allele_counting/alleleCount_info.txt"    
    """
}
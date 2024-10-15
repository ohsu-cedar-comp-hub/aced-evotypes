#!/usr/bin/env nextflow

process AlleleCounter {
    // module load R/4.3.1 //make env with R and required packages
    // module load tools/env/proxy
    container "${params.container_R}"
    publishDir "${params.bucket}/${params.sample_id}/battenberg", mode:'copy'

    input:
    path normal_bam
    path tumor_bam

    output:
    path "alleleCount_info.txt"

    script:
    """
    # get work dir directory for counting:
    workdir="$countDir/AC_${tumour_normal_id}/"
    mkdir -p ${workdir}

    R CMD BATCH "--no-restore-data --no-save --args \
    -t ${tumour_id} \
    -n ${normal_id} \
    --nb ${normal_bam_file} \
    --tb ${tumour_bam_file} \
    --cpu ${cpu} \
    -o ${workdir} " \
    runAC_forBattenberg.R "${workdir}/alleleCount_info.txt"    
    """
}
#!/usr/bin/env nextflow

process AlleleCounter {
    
    input:


    output:

    script:
    """
    Rscript '$baseDir/bin/runAC_forBattenberg.R'
    """
}
#!/usr/bin/env nextflow

// Nextflow Pipeline Version
params.release= "v0.1.1"
params.releasedate = "12-05-2024"
params.githublink = "https://github.com/ohsu-cedar-comp-hub/aced-evotypes/tree/v0.1.1/nextflow"

// import modules
include { CGPMAP } from '../tools/cgpmap.nf'
include { AlignmentMetrics; WgsMetrics; InsertMetrics } from '../tools/picard_stats.nf'

// workflow
workflow {
    // Create Channels
    basename_ch = Channel.of(params.file_1_basename, params.file_2_basename)
    path_ch = Channel.of(params.file_1_path, params.file_2_path)
    // realign using bwa-mem2
    CGPMAP(basename_ch, path_ch, params.core_zip, params.bwa_zip)

    // alignment metrics using Picard
    AlignmentMetrics(CGPMAP.out.bam, params.genome)
    WgsMetrics(CGPMAP.out.bam, params.genome)
    InsertMetrics(CGPMAP.out.bam)
}
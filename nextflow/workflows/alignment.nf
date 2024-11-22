#!/usr/bin/env nextflow

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
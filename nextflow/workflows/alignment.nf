#!/usr/bin/env nextflow

// Create Channel
filename_ch = Channel.of(params.normal, params.tumor)
basename_ch = Channel.of(params.normal_basename, params.tumor_basename)
path_ch = Channel.of(params.normal_path, params.tumor_path)

// import modules
include { CGPMAP } from '../tools/cgpmap.nf'
include { AlignmentMetrics; WgsMetrics; InsertMetrics } from '../tools/picard_stats.nf'

// workflow
workflow {
    // realign using bwa-mem2
    CGPMAP(filename_ch, basename_ch, path_ch)

    // alignment metrics using Picard
    AlignmentMetrics(CGPMAP.out.bam)
    WgsMetrics(CGPMAP.out.bam)
    InsertMetrics(CGPMAP.out.bam)
}
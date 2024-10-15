#!/usr/bin/env nextflow

// Create Channel
filename_ch = Channel.of(params.normal, params.tumor)
index_ch = Channel.of(params.normal_index, params.tumor_index)

// import modules
include { CGPWGS } from '../tools/cgpwgs.nf'
include { DellyCall, DellyTable, DellyFilter, DellyQualityFilter } from '../tools/delly.nf'
include { GRIDSS, GridssSomaticFilter } from '../tools/gridss.nf'
include { ConfigManta, MantaRunWorkflow, MantaQualityFilter } from '../tools/manta.nf'
include { ConfigStrelka, StrelkaRunWorkflow, StrelkaIndelFilter, StrelkaSnvFilter } from '../tools/strelka.nf'
include { GETPILEUPSUMMARIES, CALCULATECONTAMINATION, MUTECT2, MERGE, MERGESTATS, LEARNORIENTATION, FILTERMUTECT } from '../tools/mutect.nf'
include { AlleleCounter } from '../tools/battenberg.nf'

// workflow
workflow {
    // run CGPWGS
    CGPWGS(filename_ch, index_ch)

    // run Delly
    DellyCall(filename_ch)
    DellyTable(DellyCall.out.bcf)
    DellyFilter()
    DellyQualityFilter()

    // run GRIDSS
    GRIDSS(filename_ch)
    GridssSomaticFilter(GRIDSS.out.vcf)

    // run Manta
    ConfigManta(filename_ch)
    MantaRunWorkflow(ConfigManta.out.manta_results)
    MantaQualityFilter(MantaRunWorkflow.out.results)

    // run Strelka
    ConfigStrelka(filename_ch)
    StrelkaRunWorkflow()
    StrelkaIndelFilter()
    StrelkaSnvFilter()

    // run Mutect
    GETPILEUPSUMMARIES(filename_ch)
    CALCULATECONTAMINATION()
    MUTECT2()
    MERGE()
    MERGESTATS()
    LEARNORIENTATION()
    FILTERMUTECT()

    // run Battenberg
    AlleleCounter(filename_ch)

}
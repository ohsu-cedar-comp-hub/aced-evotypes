#!/usr/bin/env nextflow

// Nextflow Pipeline Version
params.release= "v0.1.0"
params.releasedate = "11-22-2024"
params.githublink = "https://github.com/ohsu-cedar-comp-hub/aced-evotypes/tree/v0.1.0/nextflow"

// import modules
include { SampleID } from '../tools/sample_id.nf'
include { CGPWGS } from '../tools/cgpwgs.nf'

// workflow
workflow {
    SampleID("${params.scratch_root}/${params.number_cgpwgs}/files")

    // run CGPWGS
    CGPWGS("${params.scratch_root}/${params.number_cgpwgs}/files", SampleID.out.normal_file, SampleID.out.tumor_file, params.core_zip, params.vagrent_zip, params.snv_zip, params.cnv_zip, params.qc_zip)
}
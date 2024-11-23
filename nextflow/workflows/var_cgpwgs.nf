#!/usr/bin/env nextflow

// import modules
include { SampleID } from '../tools/sample_id.nf'
include { CGPWGS } from '../tools/cgpwgs.nf'

// workflow
workflow {
    SampleID("${params.scratch_root}/${params.number_cgpwgs}/files")

    // run CGPWGS
    CGPWGS("${params.scratch_root}/${params.number_cgpwgs}/files", SampleID.out.normal_file, SampleID.out.tumor_file, params.core_zip, params.vagrent_zip, params.snv_zip, params.cnv_zip, params.qc_zip)
}
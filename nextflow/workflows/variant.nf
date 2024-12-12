#!/usr/bin/env nextflow

// Nextflow Pipeline Version
params.release= "v0.1.2"
params.releasedate = "12-12-2024"
params.githublink = "https://github.com/ohsu-cedar-comp-hub/aced-evotypes/tree/v0.1.2/nextflow"

// import modules
include { SampleID } from '../tools/sample_id.nf'
include { ConfigManta; MantaRunWorkflow; MantaQualityFilter } from '../tools/manta.nf'
include { ConfigStrelka; StrelkaRunWorkflow; StrelkaIndelFilter; StrelkaSnvFilter } from '../tools/strelka.nf'
include { DellyCall; DellyTable; DellyFilter; DellyQualityFilter } from '../tools/delly.nf'
include { GRIDSS; GridssSomaticFilter } from '../tools/gridss.nf'
include { GetPileupSummaries; CalculateContamination; MUTECT2; MERGE; MergeStats; LearnOrientation; FilterMutect} from '../tools/mutect.nf'
include { AlleleCounter } from '../tools/battenberg.nf'

// workflow
workflow {
    SampleID("${params.scratch_root}/${params.number_var}/files")

    // run Manta
    ConfigManta("${params.scratch_root}/${params.number_var}/files", SampleID.out.normal_file, SampleID.out.tumor_file, params.core_ref)
    MantaRunWorkflow("${params.scratch_root}/${params.number_var}/files", ConfigManta.out.manta, params.core_ref)
    MantaQualityFilter(MantaRunWorkflow.out.somaticSV, SampleID.out.normal_SM, SampleID.out.tumor_SM)

    // run Strelka
    ConfigStrelka("${params.scratch_root}/${params.number_var}/files", SampleID.out.normal_file, SampleID.out.tumor_file, params.core_ref)
    StrelkaRunWorkflow("${params.scratch_root}/${params.number_var}/files", ConfigStrelka.out.strelka, params.core_ref)
    StrelkaIndelFilter(StrelkaRunWorkflow.out.somatic_indels, SampleID.out.normal_SM, SampleID.out.tumor_SM)
    StrelkaSnvFilter(StrelkaRunWorkflow.out.somatic_snvs, SampleID.out.normal_SM, SampleID.out.tumor_SM)
    
    // run Delly
    DellyCall("${params.scratch_root}/${params.number_var}/files", SampleID.out.normal_file, SampleID.out.normal_SM, SampleID.out.tumor_file, SampleID.out.tumor_SM, params.genome, params.hg38_excl)
    DellyTable(DellyCall.out.bcf)
    DellyFilter(DellyCall.out.delly, DellyTable.out.samples, SampleID.out.normal_SM, SampleID.out.tumor_SM)
    DellyQualityFilter(DellyFilter.out.Bdelly, SampleID.out.normal_SM, SampleID.out.tumor_SM)

    // run GRIDSS
    GRIDSS("${params.scratch_root}/${params.number_var}/files", SampleID.out.normal_file, SampleID.out.normal_SM, SampleID.out.tumor_file, SampleID.out.tumor_SM, params.core_ref)
    GridssSomaticFilter(GRIDSS.out.vcf, params.gridss_dir, params.pondir, SampleID.out.normal_SM, SampleID.out.tumor_SM)

    // run Mutect
    GetPileupSummaries("${params.scratch_root}/${params.number_var}/files", SampleID.out.normal_file, SampleID.out.tumor_file, params.small_exac, params.small_exac_tbi, params.genome_modified)
    CalculateContamination(GetPileupSummaries.out.normal_pileup, GetPileupSummaries.out.tumor_pileup, SampleID.out.normal_SM, SampleID.out.tumor_SM)
    // Define a channel that reads each line from the intervals file. Pass this channel instead of an intervals file for splitting Mutect2 jobs
    region_ch = Channel.fromPath(params.hg38_even).splitText().map{it.trim()}
    idx_ch = Channel.of( 0..792 )

    MUTECT2(region_ch, idx_ch, "${params.scratch_root}/${params.number_var}/files", SampleID.out.normal_file, SampleID.out.normal_SM, SampleID.out.tumor_file, SampleID.out.tumor_SM, params.genome_modified, params.gnomad, params.gnomad_tbi)
    
    vcfs_channel = MUTECT2.out.vcf.collect() // collect all vcf outputs into a channel
    MERGE(vcfs_channel, params.genome_modified, SampleID.out.normal_SM, SampleID.out.tumor_SM)

    stats_ch = MUTECT2.out.stats.collect() // collect all the stats outputs into a channel
    MergeStats(stats_ch, SampleID.out.normal_SM, SampleID.out.tumor_SM)

    f1r2_ch = MUTECT2.out.f1r2.collect() // collect all the f1r2files into a channel
    LearnOrientation(f1r2_ch, SampleID.out.normal_SM, SampleID.out.tumor_SM)

    FilterMutect(MERGE.out.merged_vcf, params.genome_modified, CalculateContamination.out.contamination, 
        CalculateContamination.out.segmentation, LearnOrientation.out.read_orientation, MergeStats.out.merged_stats,
        SampleID.out.normal_SM, SampleID.out.tumor_SM)

    // run allele counter for battenberg
    AlleleCounter("${params.scratch_root}/${params.number_var}/files", params.bin, params.g1000alleles, SampleID.out.normal_SM, SampleID.out.tumor_SM)
}
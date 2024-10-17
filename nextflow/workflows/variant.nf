#!/usr/bin/env nextflow

// Create channel
file_ch = Channel.of(params.normal_path, params.tumor_path)
bai_ch = Channel.of(params.normal_bai_path, params.tumor_bai_path)
base_ch = Channel.of(params.normal_basename, params.tumor_basename)

// Define a channel that reads each line from the intervals file. Pass this channel instead of an intervals file for splitting Mutect2 jobs
region_ch = Channel.fromPath(params.hg38_even, flat: true)
idx_ch = Channel.of( 0..792 )

// import modules
include { CGPWGS } from '../tools/cgpwgs.nf'
include { ConfigManta; MantaRunWorkflow; MantaQualityFilter } from '../tools/manta.nf'
include { ConfigStrelka; StrelkaRunWorkflow; StrelkaIndelFilter; StrelkaSnvFilter } from '../tools/strelka.nf'
include { DellyCall; DellyTable; DellyFilter; DellyQualityFilter } from '../tools/delly.nf'
include { GRIDSS; GridssSomaticFilter } from '../tools/gridss.nf'
include { GetPileupSummaries; CalculateContamination; MUTECT2; MERGE; MergeStats; LearnOrientation; FILTERMUTECT} from '../tools/mutect.nf'
include { AlleleCounter } from '../tools/battenberg.nf'

// workflow
workflow {
    // run CGPWGS
    CGPWGS(params.normal_bam_path, params.normal_bai_path, params.normal_bas_path,
        params.tumor_bam_path, params.tumor_bai_path, params.tumor_bas_path,
        params.core_zip, params.vagrent_zip, params.snv_zip, params.cnv_zip, params.qc_zip)

    // run Manta
    ConfigManta(params.normal_bam_path, params.normal_bai_path, params.tumor_bam_path, params.tumor_bai_path, params.core_ref)
    MantaRunWorkflow(ConfigManta.out.manta)
    MantaQualityFilter(MantaRunWorkflow.out.somaticSV)

    // run Strelka
    ConfigStrelka(params.normal_bam_path, params.normal_bai_path, params.tumor_bam_path, params.tumor_bai_path, params.core_ref)
    StrelkaRunWorkflow(ConfigStrelka.out.strelka)
    StrelkaIndelFilter(StrelkaRunWorkflow.out.somatic_indels)
    StrelkaSnvFilter(StrelkaRunWorkflow.out.somatic_snvs)
    
    // run Delly
    DellyCall(params.normal_bam_path, params.normal_bai_path,
        params.tumor_bam_path, params.tumor_bai_path,
        params.genome, params.hg38_excl)
    DellyTable(DellyCall.out.bcf)
    DellyFilter(DellyCall.out.bcf, DellyTable.out.samples)
    DellyQualityFilter(DellyFilter.out.Bdelly)

    // run GRIDSS
    GRIDSS(params.normal_bam_path, params.tumor_bam_path, params.genome)
    GridssSomaticFilter(GRIDSS.out.vcf, params.gridss_dir, params.pondir)

    // run Mutect
    GetPileupSummaries(path_ch, bai_ch, base_ch, params.small_exac, params.genome_modified)
    CalculateContamination(GetPileupSummaries.out.pileup)
    MUTECT2(region_ch, idx_ch, param.normal_bam_path, params.tumor_bam_path, params.genome_modified, params.gnomad )
    
    vcfs_channel = MUTECT2.out.vcf.collect() // collect all vcf outputs into a channel
    MERGE(vcfs_channel, params.genome_modified)

    stats_ch = MUTECT2.out.stats.collect() // collect all the stats outputs into a channel
    MergeStats(stats_ch)

    f1r2_ch = MUTECT2.out.f1r2.collect() // collect all the f1r2files into a channel
    LearnOrientation(f1r2_ch)

    FILTERMUTECT(MERGE.out.merged_vcf, params.genome_modified, CalculateContamination.out.contamination, 
        CalculateContamination.out.segmentation, LearnOrientation.out.read_orientation, MergeStats.out.merged_stats)

    // run allele counter for battenberg
    AlleleCounter(params.normal_bam, params.tumor_bam)
}
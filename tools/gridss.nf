#!/usr/bin/env nextflow

// Step 1
process GRIDSS {
    container "${container_gridss}"

    input:
    path normal_cram_realign
    path tumor_cram_realign

    output:
    path "${seq_id_t}_v_${seq_id_n}.gridss.vcf", emit: vcf
    path "${t_n}.gridss.assembly", emit: assembly

    script:
    """
    java -jar $GRIDSS_JAR gridss \
    -r ${params.ref}/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
    -o ${seq_id_t}_v_${seq_id_n}.gridss.vcf \
    -a ${t_n}.gridss.assembly \
    normal.cram tumor.cram
    """
}

process GridssSomaticFilter {
    container "${container_gridss}"

    input:
    path GRIDSS.out.vcf

    output:
    path "${seq_id_t}_v_${seq_id_n}.gridss.filtered.vcf", emit: filtered_vcf
    path "${seq_id_t}_v_${seq_id_n}.gridss.semifiltered.vcf", emit: output_full

    script:
    """
    {rscript} /script/dir/gridss_somatic_filter \
    --pondir gridss-2.13.2 \
    --ref ${params.ref}/BSgenome.Hsapiens.UCSC.hg38 \
    --input $input \
    --output ${seq_id_t}_v_${seq_id_n}.gridss.filtered.vcf \
    --fulloutput ${seq_id_t}_v_${seq_id_n}.gridss.semifiltered.vcf \
    --normalordinal 1 \
    --tumourordinal 2 \
    --s "/script/dir/gridss_trial_17_5" \
    --gc

    """
}
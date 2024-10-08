#!/usr/bin/env nextflow

// create a channel to run realignment for normal and tumor crams
process CGPMAP {
    container "${params.container_cgpmap}"

    input: 
    path normal_file
    path tumor_file

    output:
    path "${FILE_BASENAME}.realign.OHSU", emit: realign_normal_cram
    path "${FILE_BASENAME}.realign.OHSU", emit: realign_tumor_cram

    script:
    """
    normal_basename=${params.normal_file}

    ds-cgpmap.pl \
    -r ${params.ref}/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
    -i ${params.ref}/bwa_idx_GRCh38_hla_decoy_ebv_bwamem2.tar.gz \
    -s ${normal_basename}.realign.OHSU \
    -c \
    -bwakit \
    -bm2 \
    -t 6 \
    /var/spool/data/<$FILE>
    """
}
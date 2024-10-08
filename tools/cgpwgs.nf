#!/usr/bin/env nextflow

process CGPWGS {
    container "${params.container_cgpwgs}"

    input: 
    path normal_cram_realign
    path tumor_cram_realign

    output:
    path ${params.bucket}/${params.sample_id}/CGPWGS/results
    script:
    """
    ds-cgpwgs.pl \
    -c 32 \
    -r ${params.ref}/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
    -a ${params.ref}/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
    -si ${params.ref}/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
    -cs ${params.ref}/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz \
    -qc ${params.ref}/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz \
    -e 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV' \
    -t /var/spool/data/<TUMOR CRAM> \
    -tidx /var/spool/data/<TUMOR CRAM.CRAI> \
    -n /var/spool/data/<NORMAL CRAM> \
    -nidx /var/spool/data/<NORMAL CRAM.CRAI> \
    -o /var/spool/results

    move results to bucket
    mc cp /var/spool/results ${params.bucket}/${params.sample_id}/CGPWGS/results
    """
}
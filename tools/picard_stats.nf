#!/usr/bin/env nextflow

//collect alignment metrics
process AlignmentMetrics {
    container ${params.container_picard}

    input:
    path normal_cram_realign
    path tumor_cram_realign

    output: 
    path "${SCRATCH_PATH}/${NORMAL_BASENAME}_alignment_metrics.txt"

    script:
    """
    java -jar $PICARD_JAR CollectAlignmentSummaryMetrics \
    MAX_INSERT_SIZE=100000 \
    INPUT="${SCRATCH_PATH}/${NORMAL_BASENAME}.realign.OHSU.cram" \
    OUTPUT="${SCRATCH_PATH}/${NORMAL_BASENAME}_alignment_metrics.txt" \
    REFERENCE_SEQUENCE="${params.ref}/core_ref_GRCh38_hla_decoy_ebv/genome.fa" \
    VALIDATION_STRINGENCY=LENIENT \
    IS_BISULFITE_SEQUENCED=false \
    ASSUME_SORTED=true \
    STOP_AFTER=0 \
    VERBOSITY=INFO \
    QUIET=false \
    COMPRESSION_LEVEL=5 \
    MAX_RECORDS_IN_RAM=500000 \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false
    """
}

//collect WGS metrics
process WgsMetrics {
    container ${params.container_picard}

    input:
    path normal_cram_realign
    path tumor_cram_realign

    output: 
    path "${SCRATCH_PATH}/${NORMAL_BASENAME}_coverage_metrics.txt"

    script:
    """
    java -jar $PICARD_JAR CollectWgsMetrics \
    INPUT="${SCRATCH_PATH}/${NORMAL_BASENAME}.realign.OHSU.cram" \
    OUTPUT="${SCRATCH_PATH}/${NORMAL_BASENAME}_coverage_metrics.txt" \
    REFERENCE_SEQUENCE="${params.ref}/core_ref_GRCh38_hla_decoy_ebv/genome.fa" \
    VALIDATION_STRINGENCY=SILENT \
    MINIMUM_MAPPING_QUALITY=20 \
    MINIMUM_BASE_QUALITY=20 \
    COVERAGE_CAP=250 \
    STOP_AFTER=-1 \
    INCLUDE_BQ_HISTOGRAM=false \
    VERBOSITY=INFO \
    QUIET=false \
    COMPRESSION_LEVEL=5 \
    MAX_RECORDS_IN_RAM=500000 \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false
    """
}

//collect insert size metrics
process InsertMetrics {
    container ${params.container_picard}

    input:
    path normal_cram_realign
    path tumor_cram_realign

    output: 
    path "${SCRATCH_PATH}/${NORMAL_BASENAME}_insert_metrics.pdf", emit: histogram
    path "${SCRATCH_PATH}/${NORMAL_BASENAME}_insert_metrics.txt", emit: txt

    script:
    """
    java -jar $PICARD_JAR CollectInsertSizeMetrics \
    HISTOGRAM_FILE="${SCRATCH_PATH}/${NORMAL_BASENAME}_insert_metrics.pdf" \
    INPUT="${SCRATCH_PATH}/${NORMAL_BASENAME}.realign.OHSU.cram" \
    OUTPUT="${SCRATCH_PATH}/${NORMAL_BASENAME}_insert_metrics.txt" \
    VALIDATION_STRINGENCY=LENIENT \
    DEVIATIONS=10.0 \
    MINIMUM_PCT=0.05 \
    INCLUDE_DUPLICATES=false \
    ASSUME_SORTED=true \
    STOP_AFTER=0 \
    VERBOSITY=INFO \
    QUIET=false \
    COMPRESSION_LEVEL=5 \
    MAX_RECORDS_IN_RAM=500000 \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false
    """
}
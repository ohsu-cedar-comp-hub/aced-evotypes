#!/usr/bin/env nextflow

//collect alignment metrics
process AlignmentMetrics {
    // process doesn't stop on an error condition, it just reports a message notifying you of the error event
    errorStrategy 'ignore'
    publishDir "/home/groups/EllrottLab/evotypes/nextflow/tools/${params.sample_id}/bam/alignment_metrics", mode: 'copy'

    input:
    path bam
    path genome

    output: 
    path "*alignment_metrics.output"

    script:
    """
    basename=\$(basename ${bam} .bam)

    java -jar ${params.container_picard} CollectAlignmentSummaryMetrics \
    MAX_INSERT_SIZE=100000 \
    INPUT="${bam}" \
    OUTPUT="\${basename}_alignment_metrics.txt" \
    REFERENCE_SEQUENCE="${genome}" \
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
    errorStrategy 'ignore'
    publishDir "${params.bucket}/${params.sample_id}/bam/alignment_metrics", mode:'copy'

    input:
    val bam
    path genome

    output: 
    path "*coverage_metrics.output"

    script:
    """
    java -jar ${params.container_picard} CollectWgsMetrics \
    INPUT="${bam}" \
    OUTPUT="\${basename}_coverage_metrics.txt" \
    REFERENCE_SEQUENCE="${genome}" \
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
    errorStrategy 'ignore'
    publishDir "${params.bucket}/${params.sample_id}/bam/alignment_metrics", mode:'copy'

    input:
    val bam

    output: 
    path "*insert_metrics.output"

    script:
    """
    java -jar ${params.container_picard} CollectInsertSizeMetrics \
    HISTOGRAM_FILE="\${basename}_insert_metrics.pdf" \
    INPUT="${bam}" \
    OUTPUT="\${basename}_insert_metrics.txt" \
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
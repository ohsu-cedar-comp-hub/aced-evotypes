#!/bin/bash

#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=4GB
#SBATCH --time=8:00:00

# Setting Paths and Variables
SIF='/path/to/directory/that/contains/singularity/image'
RESULTS='/path/to/output/directory'
REF='/path/to/reference/directory'
DATA='/path/to/input/BAMs'

SPECIES='Human'
# Reference assembly
ASSEMBLY='GRCh38d1'
# Sequencing protocol (e.g. WGS, WXS)
PROTOCOL='WGS'

# load module
module load /etc/modulefiles/singularity/current

# Number of CPUs
CPUS=8

# Setting Paths and Variables
SIF='/path/to/directory/that/contains/singularity/image'
DATA='/path/to/input/BAMs'
REF='/path/to/reference/directory'
MANTA_ANALYSIS_PATH='/path/to/manta/output/directory'
STRELKA_ANALYSIS_PATH='path/to/strelka/output/directory'

# run Manta
/usr/bin/time -v singularity exec \
    --cleanenv \
    --workdir $MANTA_ANALYSIS_PATH \
    --home $MANTA_ANALYSIS_PATH:/home \
    --bind $REF:/ref:ro \
    --bind $DATA:/data:ro \
    --bind $MANTA_ANALYSIS_PATH:/manta_results \
    $SIF/strelka2-manta_latest.sif \
    configManta.py \
    --normalBam /data/TCGA-X8-AAAR-10A-01D-A88X-36.WholeGenome.RP-1657.bam \
    --tumorBam /data/TCGA-X8-AAAR-01A-11D-A88X-36.WholeGenome.RP-1657.bam  \
    --referenceFasta /ref/GRCh38.d1.vd1.fa \
    --runDir /manta_results

/usr/bin/time -v singularity exec \
    --cleanenv \
    --workdir $MANTA_ANALYSIS_PATH \
    --home $MANTA_ANALYSIS_PATH:/home \
    --bind $REF:/ref:ro \
    --bind $DATA:/data:ro \
    --bind $MANTA_ANALYSIS_PATH:/manta_results \
    $SIF/strelka2-manta_latest.sif \
    /manta_results/runWorkflow.py -j $CPUS

# run Strelka
/usr/bin/time -v singularity exec \
    --cleanenv \
    --workdir $STRELKA_ANALYSIS_PATH \
    --home $STRELKA_ANALYSIS_PATH:/home \
    --bind $REF:/ref:ro \
    --bind $DATA:/data:ro \
    --bind $MANTA_ANALYSIS_PATH:/manta_results:ro \
    --bind $STRELKA_ANALYSIS_PATH:/strelka_results \
    $SIF/strelka2-manta_latest.sif \
    configureStrelkaSomaticWorkflow.py \
    --normalBam /data/TCGA-X8-AAAR-10A-01D-A88X-36.WholeGenome.RP-1657.bam \
    --tumorBam /data/TCGA-X8-AAAR-01A-11D-A88X-36.WholeGenome.RP-1657.bam \
    --referenceFasta /ref/GRCh38.d1.vd1.fa \
    --indelCandidates /manta_results/results/variants/candidateSmallIndels.vcf.gz \
    --runDir /strelka_results

/usr/bin/time -v singularity exec \
    --cleanenv \
    --workdir $STRELKA_ANALYSIS_PATH \
    --home $STRELKA_ANALYSIS_PATH:/home \
    --bind $REF:/ref:ro \
    --bind $DATA:/data:ro \
    --bind $STRELKA_ANALYSIS_PATH:/strelka_results \
    $SIF/strelka2-manta_latest.sif \
    /strelka_results/runWorkflow.py -m local -j $CPUS
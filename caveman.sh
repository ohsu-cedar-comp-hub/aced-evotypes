#!/bin/bash

#SBATCH -p exacloud
#SBATCH -c 32
#SBATCH --mem-per-cpu=4GB
#SBATCH -N 1
#SBATCH --qos=very_long_jobs
#SBATCH --time=11-0

# Load module
module load /etc/modulefiles/singularity/current
# Setting CGPWGS version
export CGPWGS_VER=2.1.0 
#singularity pull docker://quay.io/wtsicgp/dockstore-cgpwgs:$CGPWGS_VER

# Setting Paths and Variables
SIF='/path/to/directory/that/contains/singularity/image'
RESULTS='/path/to/output/directory'
REF='/path/to/reference/directory'
DATA='/path/to/input/BAMs'

SPECIES='Human'
# Reference assembly
ASSEMBLY='GRCh38d1'
# Sequencing type (e.g. pulldown|exome|genome|genomic|followup|targeted|rna_seq)
SequencingType='genome'

# run CAVEMAN
/usr/bin/time -v singularity exec \
 --cleanenv \
 --workdir $RESULTS  \
 --home $RESULTS:/home  \
 --bind $REF:/var/spool/ref:ro  \
 --bind $DATA:/var/spool/data:ro  \
 --bind $RESULTS:/var/spool/results  \
 $SIF/dockstore-cgpwgs_${CGPWGS_VER}.sif \
  caveman.pl \
  -outdir /var/spool/results \
  -reference /var/spool/ref/GRCh38.d1.vd1.fa.fai \
  -tumour-bam /var/spool/data/TCGA-X8-AAAR-01A-11D-A88X-36.WholeGenome.RP-1657.bam \
  -normal-bam /var/spool/data/TCGA-X8-AAAR-10A-01D-A88X-36.WholeGenome.RP-1657.bam \
  -ignore-file /var/spool/ref/HiDepth.bed \
  -tumour-cn /var/spool/ref/TCGA-X8-AAAR-01A-11D-A88X-36.cn.bed \
  -normal-cn /var/spool/ref/TCGA-X8-AAAR-10A-01D-A88X-36.cn.bed\
  -normal-contamination 0.3 \
  -species $SPECIES \
  -species-assembly $ASSEMBLY \
  -germline-indel /var/spool/ref/tumor.chr20_vs_normal.chr20_BL.germline.bed \
  -no-flagging \
  -seqType $SequencingType \
  -threads 32 \
  -limit 32

# tumour-cn, normal-cn, and normal-contamination were generated using ASCAT's output
# tumor.chr20_vs_normal.chr20_BL.germline.bed is a blank bed file. this causes the germline indel filter to not be applied during flagging.
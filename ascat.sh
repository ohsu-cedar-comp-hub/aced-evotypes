#!/bin/bash

#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --mem=20GB
#SBATCH --time=24:00:00

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
# Sequencing protocol (e.g. WGS, WXS)
PROTOCOL='WGS'

# run ASCAT
/usr/bin/time -v singularity exec \
 --cleanenv \
 --workdir $RESULTS  \
 --home $RESULTS:/home  \
 --bind $REF:/var/spool/ref:ro  \
 --bind $DATA:/var/spool/data:ro  \
 --bind $RESULTS:/var/spool/results  \
 $SIF/dockstore-cgpwgs_${CGPWGS_VER}.sif \
 ascat.pl \
 -outdir /var/spool/results \
 -tumour /var/spool/data/TCGA-X8-AAAR-01A-11D-A88X-36.WholeGenome.RP-1657.bam \
 -normal /var/spool/data/TCGA-X8-AAAR-10A-01D-A88X-36.WholeGenome.RP-1657.bam \
 -reference /var/spool/ref/GRCh38.d1.vd1.fa \
 -snp_gc /var/spool/ref/SnpGcCorrections.tsv \
 -protocol $PROTOCOL \
 -gender L \
 -locus /var/spool/ref/qcGenotype_GRCh38_hla_decoy_ebv/gender.tsv \
 -species $SPECIES \
 -assembly $ASSEMBLY \
 -cpus 4
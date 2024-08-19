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
export CGPWGS_VER=2.1.0 #SET ME
#singularity pull docker://quay.io/wtsicgp/dockstore-cgpwgs:$CGPWGS_VER

# cgpwgs singularity
SIF='/path/to/directory/that/contains/singularity/image'
# Path to reference directory
REF='/path/to/reference/directory'
# Path to input data directory (tumor and normal BAMs)
DATA='/path/to/input/BAMs'
# Path to output directory
RESULTS='/path/to/output/directory'

#DOWNLOADED REF FILES FROM: https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/

# run ds-cgpwgs.pl
/usr/bin/time -v singularity exec \
 --cleanenv \
 --workdir $RESULTS \
 --home $RESULTS:/home  \
 --bind $REF:/var/spool/ref:ro  \
 --bind $DATA:/var/spool/data:ro  \
 --bind $RESULTS:/var/spool/results  \
 $SIF/dockstore-cgpwgs_${CGPWGS_VER}.sif \
  ds-cgpwgs.pl \
-c 32 \
-r /var/spool/ref/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
-a /var/spool/ref/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
-si /var/spool/ref/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
-cs /var/spool/ref/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz \
-qc /var/spool/ref/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz \
-e 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV' \
-t /var/spool/data/TCGA-X8-AAAR-01A-11D-A88X-36.WholeGenome.RP-1657.bam \
-tidx /var/spool/data/TCGA-X8-AAAR-01A-11D-A88X-36.WholeGenome.RP-1657.bai \
-n /var/spool/data/TCGA-X8-AAAR-10A-01D-A88X-36.WholeGenome.RP-1657.bam \
-nidx /var/spool/data/TCGA-X8-AAAR-10A-01D-A88X-36.WholeGenome.RP-1657.bai \
-o /var/spool/results
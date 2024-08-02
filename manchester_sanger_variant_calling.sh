#!/bin/bash --login
#$ -cwd
#$ -pe smp.pe 16
#$ -t 111-116
#$ -N cgpwgs_NOORANI_corrected_111-116
set -x
set -e

#Load module - get version of apptainer that csf3 recommends
module load apps/gcc/apptainer/1.0.3

# For timings
source /mnt/bmh01-rds/UoOxford_David_W/shared/Data_for_Jamie/Script/Functions_Library.sh

# cgpwgs singularity
cgpwgs='/mnt/bmh01-rds/UoOxford_David_W/shared/cgpwgs210_singularity/dockstore-cgpwgs_2.1.0.sif'

# paths to references
refspath=/mnt/bmh01-rds/UoOxford_David_W/shared/cgpwgsHg38Refs


# Have everything for NOORANI in one output folder
AlignedGenomeDir="/mnt/bmh01-rds/UoOxford_David_W/s99384ml/oesophageal_BAM_FILES/alignment_output_Noorani"
# Get the ids of bams that Ginny says match
mapfile -t tumor_normal < "/mnt/bmh01-rds/UoOxford_David_W/s99384ml/oesophageal_BAM_FILES/CGPWGS_Jobs_28_3_24/noorani_t_vs_n.txt"
INDEX=$((SGE_TASK_ID-1))

start=$(date +%s)

# get info from the NOORANI tumour vs normal list

tumorx=$(echo ${tumor_normal[$INDEX]}| awk -F "_vs_" '{print $1}')

normalx=$(echo ${tumor_normal[$INDEX]}| awk -F "_vs_" '{print $2}')

# want the patient id to be the normal - assuming only one normal
PatientID=${normalx}

# paths to aligned BAMs
nmdatapath=$AlignedGenomeDir/$normalx
tmdatapath=$AlignedGenomeDir/$tumorx

# the workplace directory
workplace_base="/mnt/bmh01-rds/UoOxford_David_W/s99384ml/oesophageal_BAM_FILES/redone_CGPWGS_NOORANI/"
workplace=$workplace_base$PatientID/tumour_$tumorx
mkdir -p $workplace

# main job 
export  APPTAINERENV_CPU=$NSLOTS

echo "Patient=$PatientID, Normal=$normalx, Cancer=$tumorx"
echo "Calling start at "`date`

apptainer exec \
--cleanenv \
--workdir ${workplace} \
--home ${workplace}:/home \
--bind ${refspath}:/var/spool/ref:ro \
--bind ${nmdatapath}:/var/spool/nmdata:ro \
--bind ${tmdatapath}:/var/spool/tmdata:ro \
--bind ${workplace}:/var/spool/results \
$cgpwgs \
ds-cgpwgs.pl \
-c $NSLOTS \
-r /var/spool/ref/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
-a /var/spool/ref/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
-si /var/spool/ref/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
-cs /var/spool/ref/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz \
-qc /var/spool/ref/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz \
-e 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV' \
-t /var/spool/tmdata/${tumorx}.bam \
-tidx /var/spool/tmdata/${tumorx}.bam.bai \
-n /var/spool/nmdata/${normalx}.bam \
-nidx /var/spool/nmdata/${normalx}.bam.bai \
-o /var/spool/results

echo "Calling end at "`date`
echo 'Variant calling ${tumor_normal[$INDEX]} completed...'
end=$(date +%s)
elapse_time `expr $end - $start` 'Variant calling ' ${tumor_normal[$INDEX]}' done'

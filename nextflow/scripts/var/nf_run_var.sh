#!/bin/bash

#SBATCH -p batch
#SBATCH --account=cedar
#SBATCH -N 1
#SBATCH -c 25
#SBATCH --mem=120GB
#SBATCH --qos=long_jobs
#SBATCH --time=10-0
#SBATCH --output=/path/to/slurm/directory/%x.%j.out ## %j=job ID; %x=job name

# variables passed from wrapper script
# $1 = file_1
# $2 = file_2
# $3 = case_id
# $4 = scratch_root
# $5 = number_var
# $6 = params

# Activate conda environment with nextflow
echo ""
echo loading conda
conda activate nextflow

# Files
workflow="/path/to/<variant.nf>"
config="/path/to/<nf_var.config>"

# define environment variables 
echo ""
echo "************************************"
NFOUTPUT="/path/to/<nextflow_logs_variant>"
echo nextflow logs output: $NFOUTPUT

# make directory to hold input files
mkdir -p $4/$5/files
# copying normal files from bucket to scratch
mc cp rgw/gdc-esca/$3/bam/$1 $4/$5/files
mc cp rgw/gdc-esca/$3/bam/${1}.bai $4/$5/files
mc cp rgw/gdc-esca/$3/bam/${1}.bas $4/$5/files
# copying tumor files from bucket to scratch
mc cp rgw/gdc-esca/$3/bam/$2 $4/$5/files
mc cp rgw/gdc-esca/$3/bam/${2}.bai $4/$5/files
mc cp rgw/gdc-esca/$3/bam/${2}.bas $4/$5/files

# move to workdir
workdir="$4/$5"
cd $workdir

echo "Files in scratch: $workdir"
tree $workdir

# Run nextflow and force exit status 0 (|| true) to continue with tarring, exporting, and deleting the work directory regardless of run completion
echo ""
echo "************************************"
echo running nextflow

srun /usr/bin/time -v \
nextflow run $workflow \
-c $config \
-params-file $6 \
-with-report \
-with-trace || true

echo "checking what files are in the current working dir"
pwd
ls -lah

echo ""
echo  "************************************"
echo "moving slurm job id ($SLURM_JOB_ID) workdir to $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
# create output directory for this job
mkdir -p $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
# copy nextflow logs 
rsync -a .nextflow.log $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
# copy HTML execution report
rsync -a *.html $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
# # copy trace file
rsync -a *.txt $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}

# copy work directory (.command*, .exitcode, .nextflow.log)
for dir in work/*/; do
    # Skip if the directory name starts with "stage"
    if [[ "$(basename "$dir")" == stage* ]]; then
        continue
    fi
    for subdir in $dir/*/; do
        # Get the directory name without the path
        dirname=$(basename "$subdir")
        # Create a sub-directory for each task
        mkdir -p $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}/$dir
        mkdir -p $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}/$dir/$dirname
        # Sync hidden files
        rsync -a --include='.*' --exclude='*' "$subdir" "$NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}/$dir/$dirname/"
        
        # Save log files for GRIDSS (contain error messages)
        for i in $subdir/*; do
            if [[ $i == "" ]]; then
                continue
            else
                rsync -a --include='gridss.*' --exclude='*' "$i" "$NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}/$dir/$dirname/"
            fi
        done
    done
done

# remove scratch
echo "deleting workdir: $workdir"
rm -rf $workdir
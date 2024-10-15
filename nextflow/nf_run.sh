#!/bin/bash

#SBATCH --job-name=<SAMPLE_ID>_<align or variant>
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 9
#SBATCH --mem=<>
#SBATCH --time=<dd-hh:mm:ss>
#SBATCH --gres disk:2048 ## Requesting 2TB of scratch space
#SBATCH --output=</path/to/slurm_alignment OR slurm_variant/%x.%j.out ## %j=job ID; %x=job name  ## dir to store slurm output 

# Activate conda environment with nextflow
echo ""
echo loading conda
conda activate nextflow

# define environment variables 
echo ""
echo "************************************"
NFOUTPUT="/path/to/nf/logs/<nextflow_logs_alignment OR nextflow_logs_variant>"
echo nextflow logs output: $NFOUTPUT

# Create working directory in scratch
echo ""
echo "************************************"
echo making working directory
srun /usr/local/bin/mkdir-scratch.sh
SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
cd $SCRATCH_PATH

# Input files
config = $(realpath < nf_align.config OR nf_variant.config >)
workflow = $(realpath < alignment.nf OR variant.nf >)

# Run nextflow and force exit status 0 (|| true) to continue with tarring, exporting, and deleting the work directory regardless of run completion
echo ""
echo "************************************"
echo running nextflow
# module load /etc/modulefiles/singularity/current
srun /usr/bin/time -v \
nextflow run $workflow \
-c $config \
-params-file $1 \
-w $SCRATCH_PATH/work || true

# Pause shell script before continuing
sleep 30s # Wait 30 seconds

# Saving nextflow logs
echo ""
echo  "************************************"
echo "moving slurm job id ($SLURM_JOB_ID) workdir to $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
# create output directory for this job
mkdir -p $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
# copy nextflow logs 
rsync -a .nextflow.log $NFOUTPUT/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
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
    done
done

# Back out of the working directory before running the rmdir-scratch.sh script. This is only recommended on single-node jobs.
cd ../
# clean up work directory from /mnt/scratch
echo ""
echo "************************************"
echo cleaning workdir
srun /usr/local/bin/rmdir-scratch.sh
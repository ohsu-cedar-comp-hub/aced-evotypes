#!/bin/bash

#SBATCH -p batch
#SBATCH --account=cedar
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1GB
#SBATCH --time=1:00:00

json_dir="/path/to/directory/with/json/files"
script="/path/to/<nf_run_align.sh>"

for params in $json_dir/*
do
    case_id=$(jq .case_id $params | tr -d '"')
    jobname=${case_id}_align_evotypes
    sbatch -J $jobname $script $params
    # Pause shell script before continuing
    sleep 1m # Waits 1 minute
done
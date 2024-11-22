#!/bin/bash

#SBATCH -p batch
#SBATCH --account=cedar
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1GB
#SBATCH --time=2:00:00

json_dir="/path/to/directory/with/json/files"
script="/path/to/<nf_run_cgpwgs.sh>"

for params in $json_dir/*
do
    file_1=$(jq .file_1 $params | tr -d '"')
    file_2=$(jq .file_2 $params | tr -d '"')
    case_id=$(jq .case_id $params | tr -d '"')
    scratch_root=$(jq .scratch_root $params | tr -d '"')
    number=$(jq .number_cgpwgs $params | tr -d '"')
    jobname=${case_id}_cgpwgs_evotypes

    sbatch -J $jobname $script $file_1 $file_2 $case_id $scratch_root $number $params

    # Pause shell script before continuing
    sleep 1m # Waits 1 minute
done
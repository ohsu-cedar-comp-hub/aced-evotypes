#!/bin/bash

#SBATCH --job-name=gen_params
#SBATCH -p batch
#SBATCH --account=cedar
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1GB
#SBATCH --time=3:00:00

# Activate python env that contains required packages
conda activate python_3.12.4

script="/path/to/script/directory"
manifest="/path/to/sample/list"
global_params="/path/to/<global_params.json>"
output="/path/to/output/directory"

/usr/bin/time -v $script/generate_params.py \
-m $manifest \
-g $global_params \
-o $output

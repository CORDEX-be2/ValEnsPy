#!/bin/bash
#PBS -l nodes=1:ppn=128
#PBS -l walltime=06:00:00
#PBS -A 2022_205
#PBS -m abe
#PBS -N CORDEX_eval

cd $run_dir

# Activate the conda environment
source activate valenspy_env

# Run the script
python master_script.py
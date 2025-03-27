#!/bin/bash

#Get the script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

qsub -v run_dir=$DIR master_pbs_job.sh

#!/bin/bash

# Basic parameters
#PBS -N radclim_to_nc           ## Job name
#PBS -l nodes=1:ppn=1     ## 1 node, 2 processors per node (ppn=all to get a full node)
#PBS -l walltime=12:00:00 ## Max time your job will run (no more than 72:00:00)

# Situational parameters: remove one '#' at the front to use
#PBS -m abe               ## Email notifications (abe=aborted, begin and end)

#PBS -o log.out        ## Output log
#PBS -e log.err        ## Error log

module load vsc-mympirun
module load wrf-python/1.3.4.1-foss-2023a
module load h5py/3.9.0-foss-2023a
export OMPI_MCA_btl='^uct,ofi'
export OMPI_MCA_pml='ucx'
export OMPI_MCA_mtl='^ofi'

cd $PBS_O_WORKDIR 

python radclim_daily.py


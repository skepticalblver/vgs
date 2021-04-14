#!/bin/bash -f
#===============================================================================
# USERDEFINED
# This is where the batch submission is set.  The above code computes
# the total number of tasks, nodes, and other things that can be useful
# here.  Use PBS, BSUB, or whatever the local environment supports.
#===============================================================================

#SBATCH -A b1068
#SBATCH -p buyin
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -t 03:59:00
##SBATCH -m abe
##SBATCH -M howard@earth.northwesternrth.edu

cd $PBS_O_WORKDIR
module purge all
module load python

./run_earth.py

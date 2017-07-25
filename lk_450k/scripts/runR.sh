#!/bin/bash 

#SBATCH -n 1 #Number of cores 

#SBATCH -t 120 #Runtime in minutes 

#SBATCH -p serial_requeue #Partition to submit to 

#SBATCH --mem-per-cpu=100000 #Memory per cpu in MB (see also --mem)
module load centos6/R-3.0.2
R --no-save <process_450k.R

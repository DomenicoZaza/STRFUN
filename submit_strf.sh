#!/bin/bash
#SBATCH --job-name=strf1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --output=strf1.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=domenico.zaza@polito.it
#SBATCH --account=IscrB_SONORA
#SBATCH --partition=boost_usr_prod
#SBATCH --time=05:00:00

module load gcc
module load fftw/3.3.10--openmpi--4.1.6--gcc--12.2.0
module load openblas/0.3.24--gcc--12.2.0

srun --cpu-bind=cores -m block:block ./STRFUN.x

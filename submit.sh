#!/bin/bash
#SBATCH --job-name=post
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --output=post.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=domenico.zaza@polito.it
#SBATCH --account=IscrB_ALPCF_1
#SBATCH --partition=boost_usr_prod
#SBATCH --time=03:00:00

module load gcc
module load fftw/3.3.10--openmpi--4.1.6--gcc--12.2.0
module load openblas/0.3.24--gcc--12.2.0

srun --cpu-bind=cores -m block:block ./POST.x

#!/bin/bash
#SBATCH --job-name="mmgw12-nbody-10"
#SBATCH -o nbody-10
#SBATCH -e nbody-10.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=fatema.alkhanaizi@durham.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -p par6.q

source /etc/profile.d/modules.sh

module load intel/xe_2018.2

./solution-step4.o 100 1e4 5e-9 0 0 0 0 0 39.948 0 0 0 1e-12 0 0 39.948 > step4.txt
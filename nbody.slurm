#!/bin/bash
#SBATCH --job-name="mmgw12-nbody-10"
#SBATCH -o nbody-10
#SBATCH -e nbody-10.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=fatema.alkhanaizi@durham.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -p par7.q

source /etc/profile.d/modules.sh

module load intel/xe_2018.2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./solution-step2.o `python2 params.py 10`

gprof ./solution-step2.o > step2-10.txt

./solution-step6.o `python2 params.py 10`

gprof ./solution-step6.o > step6-10.txt
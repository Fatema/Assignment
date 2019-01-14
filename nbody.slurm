#!/bin/bash
#SBATCH --job-name="mmgw12-nbody"
#SBATCH -o nbody
#SBATCH -e nbody.err
#SBATCH -t 00:05:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=fatema.alkhanaizi@durham.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -p par7.q

source /etc/profile.d/modules.sh

module load intel/xe_2018.2

/ddn/apps/Cluster-Apps/intel/xe_2018.2/advisor_2018.2.0.551025/bin64/advixe-cl -collect survey -project-dir /ddn/home/mmgw12/nbody -- /ddn/home/mmgw12/Assignment/solution-step2.o `python2 params.py 1000`
/ddn/apps/Cluster-Apps/intel/xe_2018.2/advisor_2018.2.0.551025/bin64/advixe-cl -collect tripcounts -flop -project-dir /ddn/home/mmgw12/nbody -- /ddn/home/mmgw12/Assignment/solution-step2.o `python2 params.py 1000`

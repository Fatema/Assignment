#!/bin/bash
#SBATCH --job-name="mmgw12-omp"
#SBATCH -o omp
#SBATCH -e omp.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=fatema.alkhanaizi@durham.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -p par7.q

source /etc/profile.d/modules.sh

module load intel/xe_2018.2

export OMP_NUM_THREADS=1

./temp-omp.o `python2 params.py 10`
./temp-omp.o `python2 params.py 100`
./temp-omp.o `python2 params.py 1000`
./temp-omp.o `python2 params.py 10000`
./temp-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=2

./temp-omp.o `python2 params.py 10`
./temp-omp.o `python2 params.py 100`
./temp-omp.o `python2 params.py 1000`
./temp-omp.o `python2 params.py 10000`
./temp-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=4

./temp-omp.o `python2 params.py 10`
./temp-omp.o `python2 params.py 100`
./temp-omp.o `python2 params.py 1000`
./temp-omp.o `python2 params.py 10000`
./temp-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=8

./temp-omp.o `python2 params.py 10`
./temp-omp.o `python2 params.py 100`
./temp-omp.o `python2 params.py 1000`
./temp-omp.o `python2 params.py 10000`
./temp-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=16

./temp-omp.o `python2 params.py 10`
./temp-omp.o `python2 params.py 100`
./temp-omp.o `python2 params.py 1000`
./temp-omp.o `python2 params.py 10000`
./temp-omp.o `python2 params.py 100000`

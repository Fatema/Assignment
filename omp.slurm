#!/bin/bash
#SBATCH --job-name="mmgw12-omp"
#SBATCH -o omp
#SBATCH -e omp.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=fatema.alkhanaizi@durham.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -p par6.q

source /etc/profile.d/modules.sh

module load intel/xe_2018.2

export OMP_NUM_THREADS=1

./test-omp.o `python2 params.py 10`
./test-omp.o `python2 params.py 100`
./test-omp.o `python2 params.py 1000`
./test-omp.o `python2 params.py 10000`
./test-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=2

./test-omp.o `python2 params.py 10`
./test-omp.o `python2 params.py 100`
./test-omp.o `python2 params.py 1000`
./test-omp.o `python2 params.py 10000`
./test-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=4

./test-omp.o `python2 params.py 10`
./test-omp.o `python2 params.py 100`
./test-omp.o `python2 params.py 1000`
./test-omp.o `python2 params.py 10000`
./test-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=8

./test-omp.o `python2 params.py 10`
./test-omp.o `python2 params.py 100`
./test-omp.o `python2 params.py 1000`
./test-omp.o `python2 params.py 10000`
./test-omp.o `python2 params.py 100000`

export OMP_NUM_THREADS=16

./test-omp.o `python2 params.py 10`
./test-omp.o `python2 params.py 100`
./test-omp.o `python2 params.py 1000`
./test-omp.o `python2 params.py 10000`
./test-omp.o `python2 params.py 100000`

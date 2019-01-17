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

export OMP_NUM_THREADS=1

./solution-step6.o `python2 params.py 10`
gprof ./solution-step2.o > step6-10-1.txt
./solution-step6.o `python2 params.py 100`
gprof ./solution-step2.o > step6-100-1.txt
./solution-step6.o `python2 params.py 1000`
gprof ./solution-step2.o > step6-1000-1.txt

export OMP_NUM_THREADS=2

./solution-step6.o `python2 params.py 10`
gprof ./solution-step2.o > step6-10-2.txt
./solution-step6.o `python2 params.py 100`
gprof ./solution-step2.o > step6-100-2.txt
./solution-step6.o `python2 params.py 1000`
gprof ./solution-step2.o > step6-1000-2.txt
./solution-step6.o `python2 params.py 10000`
gprof ./solution-step2.o > step6-10000-2.txt

export OMP_NUM_THREADS=4

./solution-step6.o `python2 params.py 10`
gprof ./solution-step2.o > step6-10-4.txt
./solution-step6.o `python2 params.py 100`
gprof ./solution-step2.o > step6-100-4.txt
./solution-step6.o `python2 params.py 1000`
gprof ./solution-step2.o > step6-1000-4.txt
./solution-step6.o `python2 params.py 10000`
gprof ./solution-step2.o > step6-10000-4.txt

export OMP_NUM_THREADS=8

./solution-step6.o `python2 params.py 10`
gprof ./solution-step2.o > step6-10-8.txt
./solution-step6.o `python2 params.py 100`
gprof ./solution-step2.o > step6-100-8.txt
./solution-step6.o `python2 params.py 1000`
gprof ./solution-step2.o > step6-1000-8.txt
./solution-step6.o `python2 params.py 10000`
gprof ./solution-step2.o > step6-10000-8.txt

export OMP_NUM_THREADS=16

./solution-step6.o `python2 params.py 10`
gprof ./solution-step2.o > step6-10-16.txt
./solution-step6.o `python2 params.py 100`
gprof ./solution-step2.o > step6-100-16.txt
./solution-step6.o `python2 params.py 1000`
gprof ./solution-step2.o > step6-1000-16.txt
./solution-step6.o `python2 params.py 10000`
gprof ./solution-step2.o > step6-10000-16.txt

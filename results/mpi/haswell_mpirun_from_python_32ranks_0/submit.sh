#!/bin/bash
#SBATCH -o /home/morelbt/github/phd_experiments/results/mpi/haswell_mpirun_from_python_32ranks_0/submit.sh.out
#SBATCH -B 2:8:1
#SBATCH -N 2
#SBATCH -n 32
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=compute_bound
#SBATCH -t 24:00:00

python /home/morelbt/github/phd_experiments/scripts/mpi/fake_mpi_program.py 10000 100000 32
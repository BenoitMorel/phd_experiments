#!/bin/bash
#SBATCH -o /home/morelbt/github/phd_experiments/results/multi-raxml/haswellrun/phyldog_example_supporttrees_0/myjob.out
#SBATCH -B 2:8:1
#SBATCH -N 2
#SBATCH -n 32
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=compute_bound
#SBATCH -t 24:00:00

mpirun -np 32 /home/morelbt/github/multi-raxml/deps/raxml-ng/scripts/../bin/raxml-ng-mpi /home/morelbt/github/phd_experiments/results/multi-raxml/haswellrun/phyldog_example_supporttrees_0/command.txt /home/morelbt/github/phd_experiments/results/multi-raxml/haswellrun/phyldog_example_supporttrees_0/phyldog_example_supporttrees_0.svg 0
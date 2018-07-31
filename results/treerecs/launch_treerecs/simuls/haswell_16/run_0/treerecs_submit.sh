#!/bin/bash
#SBATCH -o /hits/basement/sco/morel/github/phd_experiments/results/treerecs/launch_treerecs/simuls/haswell_16/run_0/logs.out
#SBATCH -B 2:8:1
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=compute_bound
#SBATCH -t 2:00:00

/hits/basement/sco/morel/github/phd_experiments/../Treerecs/build/bin/Treerecs -g /hits/basement/sco/morel/github/phd_experiments/../datasets/simuls/geneTrees.newick -s /hits/basement/sco/morel/github/phd_experiments/../datasets/simuls/speciesTree.newick -o /hits/basement/sco/morel/github/phd_experiments/results/treerecs/launch_treerecs/simuls/haswell_16/run_0/treerecs_output -a /hits/basement/sco/morel/github/phd_experiments/../datasets/simuls/alignment.txt -t all --ale-evaluation -T 7 --select-best-tree -P 16
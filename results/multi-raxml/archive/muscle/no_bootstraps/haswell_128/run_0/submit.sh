#!/bin/bash
#SBATCH -o /hits/basement/sco/morel/github/phd_experiments/results/multi-raxml/muscle/no_bootstraps/haswell_128/run_0/logs.out
#SBATCH -B 2:8:1
#SBATCH -N 8
#SBATCH -n 128
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=compute_bound
#SBATCH -t 24:00:00

python3 /hits/basement/sco/morel/github/phd_experiments/../multi-raxml/multi-raxml/multi-raxml.py -a /hits/basement/sco/morel/github/phd_experiments/../datasets/eric_tannier/vectorbase_18/MUSCLE/ -o /hits/basement/sco/morel/github/phd_experiments/results/multi-raxml/muscle/no_bootstraps/haswell_128/run_0/multiraxml_run -r /hits/basement/sco/morel/github/phd_experiments/datasets/multi-raxml/option_files/raxml_global_options.txt -b 0 -c 128
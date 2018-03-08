#!/bin/bash
#SBATCH -o /home/morelbt/github/phd_experiments/results/multi-raxml/fromfastadir/haswell_split_256/ensembl_1/submit.sh.out
#SBATCH -B 2:8:1
#SBATCH -N 16
#SBATCH -n 256
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=compute_bound
#SBATCH -t 24:00:00

python /home/morelbt/github/phd_experiments/../multi-raxml/scripts/multiraxml_fromfastadir.py --split-scheduler /home/morelbt/github/phd_experiments/../raxml-ng/bin ~/github/datasets/ensembl_8880_15/fasta_files/ /home/morelbt/github/phd_experiments/results/multi-raxml/fromfastadir/haswell_split_256/ensembl_1 /home/morelbt/github/phd_experiments/../multi-raxml/examples/raxml_options.txt 256
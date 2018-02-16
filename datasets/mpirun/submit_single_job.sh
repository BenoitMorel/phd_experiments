#!/bin/bash
#SBATCH -o single_job.out
#SBATCH -B 2:8:1
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=compute_bound
#SBATCH -t 2:00:00

raxmlrepo=~/github/raxml-ng
raxmlexec=${raxmlrepo}/bin/raxml-ng-mpi

phdrepo=~/github/phd_experiments
msa=${phdrepo}/datasets/general/ENSGT00860000133679/fasta_files/ENSGT00860000133679.fa
outputdir=${phdrepo}/results/mpirun/single_job


outputdir1=$outputdir/mpirun1

rm -r $outputdir1
mkdir -p $outputdir1
mkdir -p $outputdir2

mpirun -np 8 ${raxmlexec} --msa $msa --model GTR --prefix $outputdir1/banana --seed 1 & 
wait



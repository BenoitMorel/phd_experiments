#!/bin/bash
#SBATCH -o multiple_nodes.out
#SBATCH -B 2:8:1
#SBATCH -N 4
#SBATCH -n 64
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=compute_bound
#SBATCH -t 2:00:00

phdrepo=~/github/phd_experiments
msa=${phdrepo}/datasets/general/ENSGT00860000133679/fasta_files/ENSGT00860000133679.fa
outputdir=${phdrepo}/results/mpirun/multiple_nodes

raxmlrepo=~/github/raxml-ng
raxmlexec=${raxmlrepo}/bin/raxml-ng-mpi

outputdir1=$outputdir/mpirun1
outputdir2=$outputdir/mpirun2
outputdir3=$outputdir/mpirun3
outputdir4=$outputdir/mpirun4

rm -r $outputdir1
rm -r $outputdir2
rm -r $outputdir3
rm -r $outputdir4
mkdir -p $outputdir1
mkdir -p $outputdir2
mkdir -p $outputdir3
mkdir -p $outputdir4

mpirun -np 16 ${raxmlexec} --msa $msa --model GTR --prefix $outputdir1/banana --seed 1 & 
mpirun -np 16 ${raxmlexec} --msa $msa --model GTR --prefix $outputdir2/banana --seed 1 &
mpirun -np 16 ${raxmlexec} --msa $msa --model GTR --prefix $outputdir3/banana --seed 1 &
mpirun -np 16 ${raxmlexec} --msa $msa --model GTR --prefix $outputdir4/banana --seed 1 &
wait



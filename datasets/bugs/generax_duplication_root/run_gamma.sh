rm gamma/ -rf &&  mpirun -np 2 ~/github/GeneRax/build/bin/generax -s species_tree.newick -f families_gamma.txt -p gamma --rec-model UndatedDTL

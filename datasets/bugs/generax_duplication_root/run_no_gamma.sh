rm no_gamma/ -rf &&  mpirun -np 2 ~/github/GeneRax/build/bin/generax -s species_tree.newick -f families.txt -p no_gamma --rec-model UndatedDTL


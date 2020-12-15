import os
import sys
import fam
from ete3 import Tree

def print_dataset_info(datadir):
  species_tree = Tree(fam.get_species_tree(datadir), format = 1)
  species = len(species_tree.get_leaves())
  print("Species number: " + str(species))
  families = fam.get_families_list(datadir)
  print("Families number: " + str(len(families)))

if (__name__ == "__main__"): 
  if (len(sys.argv) < 2): 
    print("Syntax: python " + os.path.basename(__file__) + " datadir")
    exit(1)
  print_dataset_info(sys.argv[1])

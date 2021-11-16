import os
import sys
import subprocess
import shutil
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
from read_tree import read_tree
import saved_metrics
import get_dico


def estimate(datadir):
  species_tree = read_tree(fam.get_true_species_tree(datadir))
  species_number = len(species_tree.get_leaves())
  print(species_number)
  cov = 0.0
  for family in fam.get_families_list(datadir):
    d = get_dico.get_species_to_genes_family(datadir, family)
    cov += len(d)
 
  cov = cov / float(len(fam.get_families_list(datadir)))
  print(cov / species_number)


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.ex(sys.argv[1], sys.argv[2], sys.argv[3])
  estimate(sys.argv[1])


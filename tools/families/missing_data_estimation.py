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
  per_species_coverage = {}
  for leaf in species_tree.get_leaves():
    per_species_coverage[leaf.name] = 0
  species_number = len(species_tree.get_leaves())
  print(species_number)
  writer = open(fam.get_species_coverage(datadir), "w")
  cov = 0.0
  for family in fam.get_families_list(datadir):
    d = get_dico.get_species_to_genes_family(datadir, family)
    for species in d:
      per_species_coverage[species] += 1

    cov += len(d)
 
  fam_number = len(fam.get_families_list(datadir))
  for species in per_species_coverage:
    c = float(per_species_coverage[species])
    print(species + " : " + str(c / float(fam_number)))
    writer.write(species + " : " + str(c / float(fam_number)) + "\n")
  cov = cov / float(fam_number)
  print(cov / species_number)


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  estimate(sys.argv[1])


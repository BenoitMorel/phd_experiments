import os
import sys
from ete3 import Tree
sys.path.insert(0, 'tools/mappings')

import fam
import get_dico

def evaluate(datadir):
  coverage = {}
  total = 0
  species_tree = Tree(fam.get_species_tree(datadir), format = 1)
  for species in species_tree.get_leaf_names():
    coverage[species] = 0

  for family in fam.get_families_list(datadir):
    gene_to_species = get_dico.get_gene_to_species(datadir, family)
    tree = Tree(fam.get_true_tree(datadir, family), format = 1)
    leaves = tree.get_leaf_names()
    covered_species = set()
    for leaf in leaves:
      covered_species.add(gene_to_species[leaf])
    for species in covered_species:
      coverage[species] += 1
    total += 1

  missing_data_file = fam.get_missing_data_file(datadir)
  with open(missing_data_file, "w") as writer:
    for species in coverage:
      fm = 1.0 - float(coverage[species]) / float(total)
      writer.write(species + " " + str(fm) + "\n")
  print("Wrote missing gene information in " + missing_data_file)

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  evaluate(sys.argv[1])

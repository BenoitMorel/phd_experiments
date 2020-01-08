import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/mappings')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import ete3
import get_dico
import rf_distance

def get_discordance_rate(datadir):
  species_tree = ete3.Tree(fam.get_species_tree(datadir), 1)
  sum_average_rf = 0.0
  for family in fam.get_families_list(datadir):
    gene_tree = ete3.Tree(fam.get_true_tree(datadir, family), 1)
    d = get_dico.get_gene_to_species(datadir, family)
    for leaf in gene_tree:
      leaf.name = d[leaf.name]
    rf = rf_distance.ete3_rf(species_tree, gene_tree)
    relative_rf = float(rf[0]) / float(rf[1])
    sum_average_rf += relative_rf
  average_relative_rf = sum_average_rf / float(len(fam.get_families_list(datadir)))
  return average_relative_rf

if (__name__== "__main__"):
  if len(sys.argv) != 2:
    print("Syntax: python " + os.path.basename(__file__) + " datadir")
    exit(1)
  datadir = sys.argv[1]
  print("Discordance rate: " + str(get_discordance_rate(datadir)))



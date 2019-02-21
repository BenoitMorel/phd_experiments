import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/stag')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import extract_gene_trees
import extract_mapping
import rf_distance
import stag

def run_stag_from_families_and_method(dataset, method_tag):
  output_species_tree = os.path.join(dataset, "stag", "stag_" + method_tag + "SpeciesTree.newick")
  families_path = os.path.join(dataset, "families")
  stag_gene_trees = extract_gene_trees.extract(dataset, method_tag)
  mapping_file = os.path.join(dataset, "stag_mapping.txt")
  if (not os.path.isfile(mapping_file)):
    extract_mapping.extract(dataset, mapping_file)
  stag.run_stag(mapping_file, stag_gene_trees, output_species_tree)
  return output_species_tree

def run_stag_from_families(dataset):
  methods = ["raxml", "treerecs", "phyldog", "notung"]
  true_species_tree = os.path.join(dataset, "speciesTree.newick")
  for method in methods:
    try:
      output_species_tree = run_stag_from_families_and_method(dataset, method)
      rf_tuple = rf_distance.get_rf_distance_tuple(output_species_tree, true_species_tree)
      print("RF distance between " + method + " and true species tree: \t" + str(rf_tuple))
    except:
      print("failed to run stag for " + dataset + " and method " + method)
  print("") 

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python run_stag.py dataset")
    sys.exit(1)
  run_stag_from_families(sys.argv[1])


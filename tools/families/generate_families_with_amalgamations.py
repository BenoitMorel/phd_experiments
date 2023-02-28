import sys
import os
import shutil
import subprocess
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree
from ete3 import Tree
import re

def  treat_amalgamation(amalgamation, datadir, family, authorized_species):
  f = open(amalgamation)
  f.readline() # skip first comment
  line = f.readline()
  print(family) 
  tree = None
  try:
    tree = Tree(line)
  except:
    print("failed to read tree, skipping " + amalgamation)
    return
  if (len(tree.get_leaf_names()) < 4):
    return
  fam.init_family_directories(datadir, family)
  exp.mkdir(fam.get_amalgamation_dir(datadir, family))

  output_mapping_file = fam.get_mappings(datadir, family)
  output_amalgamation = fam.get_amalgamation(datadir, family, "true", "true")
  shutil.copy(amalgamation, output_amalgamation)
  species_to_genes = {}
  for gene in tree.get_leaf_names():
    species = gene.split("_")[0]
    assert(species in authorized_species)
    if (not species in species_to_genes):
      species_to_genes[species] = [gene]
    else:
      species_to_genes[species].append(gene)
  with open(output_mapping_file, "w") as writer:
    for species in species_to_genes:
      genes = species_to_genes[species]
      writer.write(species + ":" + ";".join(genes) + "\n")
  

def generate(amalgamations_dir, species_tree_path, datadir):
  fam.init_top_directories(datadir)
  
  true_species_tree = fam.get_species_tree(datadir)
  authorized_species = set()
  authorized_species = Tree(species_tree_path, format=1).get_leaf_names()
  shutil.copy(species_tree_path, true_species_tree)
  for f in os.listdir(amalgamations_dir):
    family = f.split(".")[0]
    amalgamation = os.path.join(amalgamations_dir, f)
    treat_amalgamation(amalgamation, datadir, family, authorized_species)
  fam.postprocess_datadir(datadir)
  


if (__name__ == "__main__"): 
  if (len(sys.argv) != 4): 
    print("Syntax: python " + os.path.basename(__file__) + " amalgamations_dir species_tree  output_dir")
    exit(1)
  amalgamations_dir = sys.argv[1]
  species_tree_path = sys.argv[2]
  datadir = sys.argv[3]
  generate(amalgamations_dir, species_tree_path, datadir)





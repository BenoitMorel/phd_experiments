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

def  treat_gene_tree(gene_tree, datadir, family, authorized_species):
  tree = Tree(gene_tree)
  print(gene_tree)
  if (len(tree.get_leaf_names()) < 4):
    return
  fam.init_family_directories(datadir, family)
  output_mapping_file = fam.get_mappings(datadir, family)
  output_gene_tree = fam.get_true_tree(datadir, family)
  shutil.copy(gene_tree, output_gene_tree)
  with open(output_mapping_file, "w") as writer:
    for gene in tree.get_leaf_names():
      assert(gene in authorized_species)
      writer.write(gene + ":" + gene + "\n")
  

def generate(gene_trees_dir, species_tree_path, datadir):
  fam.init_top_directories(datadir)
  
  true_species_tree = fam.get_species_tree(datadir)
  authorized_species = set()
  authorized_species = Tree(species_tree_path, format=1).get_leaf_names()
  shutil.copy(species_tree_path, true_species_tree)
  for f in os.listdir(gene_trees_dir):
    family = f.split(".")[0]
    gene_tree = os.path.join(gene_trees_dir, f)
    treat_gene_tree(gene_tree, datadir, family, authorized_species)
  fam.postprocess_datadir(datadir)
  


if (__name__ == "__main__"): 
  if (len(sys.argv) != 4): 
    print("Syntax: python " + os.path.basename(__file__) + " gene_trees_dict species_tree  output_dir")
    exit(1)
  gene_trees_dir = sys.argv[1]
  species_tree_path = sys.argv[2]
  datadir = sys.argv[3]
  generate(gene_trees_dir, species_tree_path, datadir)




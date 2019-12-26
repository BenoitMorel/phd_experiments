import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import experiments as exp
import ete3

def get_tree(ale_file):
  for line in open(ale_file):
    if (";" in line):
      return line
  return None

def get_leaves(tree_file):
  return ete3.Tree(tree_file,format=1).get_leaf_names()

def create_mapping_from_gene_tree(input_tree, output_mapping):
  leaves = ete3.Tree(input_tree, format=1).get_leaf_names()
  species_to_gene = {}
  for leaf in leaves:
    species = "_".join(leaf.split("_")[1:])
    if (not species in species_to_gene):
      species_to_gene[species] = []
    species_to_gene[species].append(leaf)
  with open(output_mapping, "w") as writer:
    for species in species_to_gene:
      genes = species_to_gene[species]
      writer.write(species + ":" + ";".join(genes) + "\n")

def export(input_trees_dir, species_tree, datadir):
  
  # init directories
  print("Starts generation")
  fam.init_top_directories(datadir)
  
  # species tree
  true_species_tree = fam.get_species_tree(datadir)
  shutil.copy(species_tree, true_species_tree)
  fam.init_top_directories(datadir)
  
  # families 
  print("Init families")
  families = []
  for f in os.listdir(input_trees_dir):
    family = f.split(".")[0]
    families.append(f.split(".")[0])
  fam.init_families_directories(datadir, families)

  # fill families
  print("Fill families")
  for family in families:
    ale_file = os.path.join(input_trees_dir, family + ".ale")
    output_tree = fam.get_true_tree(datadir, family)
    tree_string = get_tree(ale_file)
    with open(output_tree, "w") as writer:
      writer.write(tree_string)
    mapping_file = fam.get_mappings(datadir, family)
    create_mapping_from_gene_tree(tree_string, mapping_file) 
    
  print("post process")
  fam.postprocess_datadir(datadir)
      

if (__name__ == "__main__"): 
  if (len(sys.argv) != 4):
    print("syntax: input_trees_dir input_species_tree datadir")
    exit(1)

  input_trees_dir = sys.argv[1]
  species_tree = sys.argv[2]
  datadir = sys.argv[3]
  export(input_trees_dir, species_tree, datadir)



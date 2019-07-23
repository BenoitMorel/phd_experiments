import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import stag
import fam
from read_tree import read_tree


def init_gene_trees_file(datadir, subst_model, output_dir):
  gene_trees_file = os.path.join(output_dir, "gene_trees.txt")
  trees = []
  for family in fam.get_families_list(datadir):
    tree = fam.get_raxml_tree(datadir, subst_model, family)
  return gene_trees_file

def init_species_file(datadir, output_dir):
  species_file = os.path.join(output_dir, "species.txt")
  leaves = read_tree(fam.get_species_tree(datadir)).get_leaves()
  with open(species_file, "w") as writer:
    writer.write("\n".join(leaves))
  return species_file

def build_guenomu_config_file(gene_trees_file, species_tree_file, output_dir):
  config_file = os.path.join(output_dir, "config.txt")
  return config_file

def execute_guenomu(dataset_dir, config_file):
  pass

def run_guenomu(datadir, subst_model):
  output_dir = fam.get_run_dir(datadir, subst_model, "guenomu_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  output_species_tree = fam.get_species_tree(datadir, subst_model, "guenomu")
  gene_trees_file = init_gene_trees_file(datadir, subst_model, output_dir)
  species_file = init_species_file(datadir, output_dir)
  config_file = build_guenomu_config_file(gene_trees_file, species_file, output_dir)
  execute_guenomu(dataset_dir, config_file)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python run_guenomu.py datadir subst_model")
    sys.exit(1)
  run_guenomu(sys.argv[1], sys.argv[2])
  


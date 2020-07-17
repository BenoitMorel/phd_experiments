import os
import sys
import subprocess
import shutil
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/stag')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import stag
import saved_metrics
import fam


def init_gene_trees_file(datadir, gene_trees, subst_model, stag_gene_trees_dir):
  exp.mkdir(stag_gene_trees_dir) 
  print("datadir " + datadir)
  for family in fam.get_families_list(datadir):
    gene_tree = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
    dest = os.path.join(stag_gene_trees_dir, family + ".newick")
    shutil.copyfile(gene_tree, dest)

def init_mappings(datadir, output_mapping_file):
  dico = {}
  for family in fam.get_families_list(datadir):
    mapping_file = fam.get_mappings(datadir, family)
    print(family)
    for line in open(mapping_file).readlines():
      sp1 = line.split(":")
      species = sp1[0]
      if (not species in dico):
        dico[species] = []
      genes = sp1[1][:-1].split(";")
      dico[species].extend(genes)
  with open(output_mapping_file, "w") as writer:
    for species in dico:
      for gene in dico[species]:
        writer.write(gene + " " + species + "\n")
  print("wrote into " + output_mapping_file)

def run_stag(datadir, gene_trees,  subst_model):
  run_name = "stag_" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  output_species_tree = fam.get_species_tree(datadir, subst_model, run_name) 
  stag_gene_trees_dir = os.path.join(output_dir, "stag_gene_trees")
  init_gene_trees_file(datadir, gene_trees, subst_model, stag_gene_trees_dir)
  mapping_file = os.path.join(output_dir, "stag_mapping.txt")
  init_mappings(datadir, mapping_file)
  start = time.time()
  stag.run_stag(mapping_file, stag_gene_trees_dir, output_species_tree)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  true_species_tree = fam.get_species_tree(datadir)
  print("Stag output: " + output_species_tree)


if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python run_stag.py datadir gene_trees subst_model")
    sys.exit(1)
  run_stag(sys.argv[1], sys.argv[2], sys.argv[3])
  

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
from read_tree import read_trees_list
import saved_metrics
import get_dico
import species_analyze

def init_gene_trees_file(datadir, gene_trees, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  print(filepath) 
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      skip_family = False
      gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
      mapping_path = fam.get_mappings(datadir, family)
      m = get_dico.get_gene_to_species(datadir, family)
      for line in open(mapping_path).readlines():
        split = line.replace("\n", "").split(":")
        species = split[0]
        genes = split[1].split(";")
        if (len(genes) > 1):
          skip_family = True
        for gene in genes:
          m[gene] = species
      if (skip_family):
        continue
      trees = read_trees_list(gene_tree_path)
      for tree in trees:
        for leaf in tree:
          leaf.name = m[leaf.name]
        writer.write(tree.write().replace("[0123456789]e-", ""))
        writer.write("\n")
  return filepath


def exec_astrid(gene_trees_file, output_species_tree_file, mode):
  command = []
  command.append(exp.astrid_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-a")
  command.append(output_species_tree_file)
  command.append("-o")
  command.append(output_species_tree_file)
  #if (mode == "bionj"):
  #  command.append("--bionj")
  if (mode == "fastme"):
    command.append("-s")
  elif(mode == "default"):
    pass
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command)#, stderr=FNULL)

def run_astrid(datadir, method, subst_model, mode):
  name = "astrid-" + mode + "_" + method
  output_dir = fam.get_run_dir(datadir, subst_model, name + "_run")
  temp_species_tree = os.path.join(output_dir, "species_tree.newick")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, method, subst_model, output_dir)
  print("Start executing astrid")
  start = time.time()
  exec_astrid(gene_trees_file, temp_species_tree, mode)
  output_species_tree = fam.get_species_tree(datadir, subst_model, name)
  shutil.move(temp_species_tree, output_species_tree)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(name, subst_model), time1, "runtimes") 

if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python run_astrid.py datadir gene_trees subst_model [mode=default,,fastme]")
    sys.exit(1)
  datadir = sys.argv[1]
  method = sys.argv[2]
  subst_model = sys.argv[3]
  mode = "fastme"
  if (len(sys.argv) == 5):
    mode = sys.argv[4]
  run_astrid(datadir, method, subst_model, mode)
  species_analyze.analyze(datadir) 
  





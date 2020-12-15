import shutil
import sys
import os
import subprocess
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/mappings')
import get_dico
from read_tree import read_trees_list
import experiments as exp
import fam
import species_analyze
import saved_metrics

def init_gene_trees_file(datadir, gene_trees, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
      mapping_path = fam.get_mappings(datadir, family)
      m = get_dico.get_gene_to_species(datadir, family)
      for line in open(mapping_path).readlines():
        split = line.replace("\n", "").split(":")
        species = split[0]
        genes = split[1].split(";")
        for gene in genes:
          m[gene] = species
      trees = read_trees_list(gene_tree_path)
      for tree in trees:
        for leaf in tree:
          leaf.name = m[leaf.name]
        writer.write(tree.write().replace("e-", ""))
        writer.write("\n")
  return filepath


def exec_paul(algo, gene_trees_file, output_tree):
  command = []
  command.append(exp.paul_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-o")
  command.append(output_tree)
  if (algo == "NJst"):
    pass
  elif (algo == "MiniNJ"):
    command.append("-m")
  elif (algo == "WMiniNJ"):
    command.append("-m")
    command.append("-w")
  elif (algo == "TagNJ"):
    command.append("-t")
  elif (algo == "RootTagNJ"):
    command.append("-t")
    command.append("-r")
  else:
    print("Invalid algorithm in exec_paul: " + algo)
    assert(False)
  print(" ".join(command))
  subprocess.check_call(command)
    

def run_paul(datadir, algo, gene_trees, subst_model):
  run_name = "paul-" + algo + "_" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  output_tree = fam.get_species_tree(datadir, subst_model, run_name) 
  start = time.time()
  exec_paul(algo, gene_trees_file, output_tree)
  time1 = (time.time() - start)
  print("Runtime: " + str(time1) + "s")
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 

  


if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_paul.py datadir algo gene_trees subst_models.")
    sys.exit(0)


  datadir = sys.argv[1]
  algo = sys.argv[2]
  gene_trees = sys.argv[3]
  subst_model = sys.argv[4] 
  run_paul(datadir, algo, gene_trees, subst_model)
  species_analyze.analyze(datadir)




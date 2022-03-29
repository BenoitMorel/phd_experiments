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
from read_tree import read_trees_list
import saved_metrics
import get_dico
import species_analyze

def init_gene_trees_file(datadir, gene_trees, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  print(filepath) 
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
        writer.write(tree.write().replace("e-", "").replace("e+",""))
        writer.write("\n")
  return filepath


def exec_stells(gene_trees_file, output_dir, cores):
  command = []
  command.append(exp.stells_exec)
  command.append("-g")
  command.append(gene_trees_file)
  if (cores > 1):
    command.append("-t")
    command.append(str(cores))
  out = os.path.join(output_dir, "out.txt")
  print(" ".join(command))
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command)
  print(res)
  lines = res.split("\n")
  key = "The newick format of the inferred MLE species tree: "
  for line in lines:
    if (line.startswith(key)):
      tree = line.replace(key, "") + ";"
      with open(out, "w") as writer:
        writer.write(tree)
  return out

def run_stells(datadir, gene_trees, subst_model, cores = 1):
  run_name = "stells_" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  start = time.time()
  stells_tree = exec_stells(gene_trees_file, output_dir, cores)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  shutil.copy(stells_tree, fam.get_species_tree(datadir, subst_model, run_name))

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model")
    sys.exit(1)
  datadir = sys.argv[1]
  gene_trees = sys.argv[2]
  subst_model = sys.argv[3]
  cores = int(sys.argv[4])
  run_stells(datadir, gene_trees, subst_model, cores)
  species_analyze.analyze(sys.argv[1])
  




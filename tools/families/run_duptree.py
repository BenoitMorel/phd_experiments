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


def exec_duptree(gene_trees_file, output_dir):
  command = []
  command.append(exp.duptree_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-o")
  command.append(os.path.join(output_dir, "plop.txt"))
  #command.append("-r")
  #command.append("all")
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command, stderr=FNULL)
  return res.split("\n")[-3]

def run_duptree(datadir, gene_trees, subst_model):
  run_name = "duptree-" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  start = time.time()
  species_tree_str = exec_duptree(gene_trees_file, output_dir)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  assert(species_tree_str.startswith("("))
  with open(fam.get_species_tree(datadir, subst_model, run_name), "w") as writer:
    writer.write(species_tree_str)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python run_duptree.py datadir gene_trees subst_model")
    sys.exit(1)
  run_duptree(sys.argv[1], sys.argv[2], sys.argv[3])
  



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
import saved_metrics
import get_dico

def init_gene_trees_file(datadir, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
 
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      raxml_tree_path = fam.get_raxml_tree(datadir, subst_model, family)
      mapping_path = fam.get_mappings(datadir, family)
      tree = read_tree(raxml_tree_path)
      m = get_dico.get_gene_to_species(datadir, family)
      for line in open(mapping_path).readlines():
        split = line.replace("\n", "").split(":")
        species = split[0]
        genes = split[1].split(";")
        for gene in genes:
          m[gene] = species
      for leaf in tree:
        leaf.name = m[leaf.name]
      writer.write(tree.write().replace("e-", ""))
      writer.write("\n")
  return filepath


def exec_fastrfs(gene_trees_file):
  command = []
  command.append(exp.fastrfs_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-o")
  command.append("plop.txt")
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command, stderr=FNULL)
  return res.split("\n")[-3]

def run_fastrfs(datadir, subst_model):
  output_dir = fam.get_run_dir(datadir, subst_model, "fastrfs_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, subst_model, output_dir)
  start = time.time()
  species_tree_str = exec_fastrfs(gene_trees_file)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name("fastrfs", subst_model), time1, "runtimes") 
  assert(species_tree_str.startswith("("))
  with open(fam.get_species_tree(datadir, subst_model, "fastrfs"), "w") as writer:
    writer.write(species_tree_str)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python run_fastrfs.py datadir subst_model")
    sys.exit(1)
  run_fastrfs(sys.argv[1], sys.argv[2])
  



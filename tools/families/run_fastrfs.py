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
from read_tree import read_trees_list
import run_astral_multi as astral

def init_gene_trees_file(datadir, gene_trees, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      skip_family = False
      gene_tree_path = fam.get_gene_tree(datadir, subst_model, family, gene_trees)
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
        writer.write(tree.write().replace("e-", ""))
        writer.write("\n")
  return filepath



def exec_fastrfs(gene_trees_file, output_prefix):
  command = []
  command.append(exp.fastrfs_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-o")
  command.append(output_prefix)
  command.append("--nostrict")
  command.append("--nomajority")
  FNULL = open(os.devnull, 'w')
  subprocess.check_call(command, stdout=FNULL, stderr=FNULL)

def extract_species_trees(datadir, gene_trees, subst_model, output_prefix, method):
  greedy = output_prefix + "." + method
  greedy_out = fam.get_species_tree(datadir, subst_model, "fastrfs-" + gene_trees + "_" + method)
  shutil.copyfile(greedy, greedy_out)


def run_fastrfs(datadir, gene_trees, subst_model):
  run_name = "fastrfs-" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  prefix = os.path.join(output_dir, "output")
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  start = time.time()
  exec_fastrfs(gene_trees_file, prefix)
  time1 = (time.time() - start)
  for method in ["greedy", "single"]:
    saved_metrics.save_metrics(datadir, fam.get_run_name("fastrfs_" + method, subst_model), time1, "runtimes") 
    extract_species_trees(datadir, gene_trees, subst_model, prefix, method)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python run_fastrfs.py datadir gene_trees subst_model")
    sys.exit(1)
  run_fastrfs(sys.argv[1], sys.argv[2], sys.argv[3])
  



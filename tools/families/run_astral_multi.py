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
      writer.write(open(raxml_tree_path).read().replace("\n", "") + "\n")
  return filepath

def init_mapping_file(datadir, output_dir):
  res = os.path.join(output_dir, "astral_mappings.txt")
  with open(res, "w") as writer:
    dico = get_dico.get_species_to_genes(datadir)
    for species in dico:
      writer.write(species)
      writer.write(":")
      writer.write(",".join(dico[species]))
      writer.write("\n")
  return res

def exec_astral(gene_trees_file, mapping_file, output_species_tree_file):
  command = []
  command.append("java")
  command.append("-jar")
  command.append(exp.astral_jar)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-a")
  command.append(mapping_file)
  command.append("-o")
  command.append(output_species_tree_file)
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command, stderr=FNULL)
  print(res)

def run_astral(datadir, subst_model):
  output_dir = fam.get_run_dir(datadir, subst_model, "astral_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, subst_model, output_dir)
  mapping_file = init_mapping_file(datadir, output_dir)
  start = time.time()
  exec_astral(gene_trees_file, mapping_file, fam.get_species_tree(datadir, subst_model, "astral"))
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name("astral", subst_model), time1, "runtimes") 

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python run_astral.py datadir subst_model")
    sys.exit(1)
  run_astral(sys.argv[1], sys.argv[2])
  




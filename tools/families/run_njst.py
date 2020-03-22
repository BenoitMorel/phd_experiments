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

def init_gene_trees_file(datadir, method, subst_model, output_dir):
  filepath = os.path.join(output_dir, "gene_trees.txt")
  with open(filepath, "w") as writer:
    for family in fam.get_families_list(datadir):
      raxml_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, method)
      writer.write(open(raxml_tree_path).read().replace("\n", "") + "\n")
  return filepath

def init_mapping_file(datadir, output_dir):
  res = os.path.join(output_dir, "njst_mappings.txt")
  with open(res, "w") as writer:
    dico = get_dico.get_species_to_genes(datadir)
    for species in dico:
      for gene in dico[species]:
        writer.write(gene + " " + species + "\n")
  return res

def exec_njst(gene_trees_file, mapping_file, output_species_tree_file, algo_type):
  command = []
  command.append("Rscript")
  command.append(exp.njstm_script)
  command.append(gene_trees_file)
  command.append(mapping_file)
  command.append(algo_type)
  command.append(output_species_tree_file)
  FNULL = open(os.devnull, 'w')
  res = subprocess.check_output(command, stderr=FNULL)
  print(res)

def run_njst_all(datadir, method, subst_model):
  run_njst(datadir, method, subst_model, "reweighted")
  run_njst(datadir, method, subst_model, "liu")
  run_njst(datadir, method, subst_model, "original")

def run_njst(datadir, method, subst_model, algo_type):
  algo_name = "njst-" + algo_type
  output_dir = fam.get_run_dir(datadir, subst_model, algo_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, method, subst_model, output_dir)
  mapping_file = init_mapping_file(datadir, output_dir)
  print("Start executing " + algo_name)
  start = time.time()
  exec_njst(gene_trees_file, mapping_file, fam.get_species_tree(datadir, subst_model, algo_name), algo_type)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(algo_name, subst_model), time1, "runtimes") 

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python run_njst.py datadir gene_trees subst_model")
    sys.exit(1)
  run_njst_all(sys.argv[1], sys.argv[2], sys.argv[3])
  




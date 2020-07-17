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
import run_astral_multi as astral
import run_duptree
import species_analyze

def init_gene_trees_file(datadir, gene_trees, subst_model, output_dir):
  
  return run_duptree.init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  
  return astral_gene_trees
  return
  print("astral_gene_trees: " + astral_gene_trees)
  astral_mappings = astral.init_mapping_file(datadir, output_dir)
  command = []
  command.append("python")
  command.append(exp.prepare_fastrfs_script)
  command.append("-i")
  command.append(astral_gene_trees)
  command.append("-a")
  command.append(astral_mappings)
  subprocess.check_call(command)
  res = astral_gene_trees[:-4] + "_g_trees-mult.trees"
  return res


def exec_fastmulrfs(gene_trees_file, output_prefix):
  preprocessed_gene_trees = gene_trees_file + ".preprocessed"
  pre_command = []
  pre_command.append("python")
  pre_command.append(exp.fastmulrfs_preprocess)
  pre_command.append("-i")
  pre_command.append(gene_trees_file)
  pre_command.append("-o")
  pre_command.append(preprocessed_gene_trees)
  subprocess.check_call(pre_command)
  
  
  command = []
  command.append(exp.fastrfs_exec)
  command.append("-i")
  command.append(preprocessed_gene_trees)
  command.append("-o")
  command.append(output_prefix)
  FNULL = open(os.devnull, 'w')
  subprocess.check_call(command)
  #subprocess.check_call(command, stdout=FNULL, stderr=FNULL)

def extract_species_trees(datadir, subst_model, output_prefix, method, sub_run_name):
  greedy = output_prefix + "." + method
  greedy_out = fam.get_species_tree(datadir, subst_model, sub_run_name)
  shutil.copyfile(greedy, greedy_out)
  print(greedy_out)

def run_fastmulrfs(datadir, gene_trees, subst_model):
  run_name = "fastmulrfs_" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  prefix = os.path.join(output_dir, "output")
  gene_trees_file = init_gene_trees_file(datadir, gene_trees, subst_model, output_dir)
  start = time.time()
  exec_fastmulrfs(gene_trees_file, prefix)
  time1 = (time.time() - start)
  for method in ["greedy", "single", "strict", "majority"]:
    sub_run_name = run_name.replace("fastmulrfs", "fastmulrfs-" + method)
    saved_metrics.save_metrics(datadir, fam.get_run_name(sub_run_name, subst_model), time1, "runtimes") 
    extract_species_trees(datadir, subst_model, prefix, method, sub_run_name)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python run_fastmulrfs.py datadir gene_trees subst_model")
    sys.exit(1)
  run_fastmulrfs(sys.argv[1], sys.argv[2], sys.argv[3])
  species_analyze.analyze(sys.argv[1])
  



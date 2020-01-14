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

def init_gene_trees_file(datadir, subst_model, output_dir):
  astral_gene_trees = astral.init_gene_trees_file(datadir, "raxml-ng", subst_model, output_dir)
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
  res = astral_gene_trees[:-4] + "-for-fastrfs.txt"
  return res


def exec_fastrfs(gene_trees_file, output_prefix):
  command = []
  command.append(exp.fastrfs_exec)
  command.append("-i")
  command.append(gene_trees_file)
  command.append("-o")
  command.append(output_prefix)
  FNULL = open(os.devnull, 'w')
  subprocess.check_call(command, stdout=FNULL, stderr=FNULL)

def extract_aux(datadir, subst_model, output_prefix, method):
  greedy = output_prefix + "." + method
  greedy_out = fam.get_species_tree(datadir, subst_model, "fastrfs_" + method)
  shutil.copyfile(greedy, greedy_out)

def extract_species_trees(datadir, subst_model, output_prefix):
  extract_aux(datadir, subst_model, output_prefix, "greedy")
  extract_aux(datadir, subst_model, output_prefix, "single")
  extract_aux(datadir, subst_model, output_prefix, "strict")
  extract_aux(datadir, subst_model, output_prefix, "majority")

def run_fastrfs(datadir, subst_model):
  output_dir = fam.get_run_dir(datadir, subst_model, "fastrfs_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  prefix = os.path.join(output_dir, "output")
  gene_trees_file = init_gene_trees_file(datadir, subst_model, output_dir)
  start = time.time()
  exec_fastrfs(gene_trees_file, prefix)
  time1 = (time.time() - start)
  for method in ["greedy", "single", "strict", "majority"]:
    saved_metrics.save_metrics(datadir, fam.get_run_name("fastrfs_" + method, subst_model), time1, "runtimes") 
  extract_species_trees(datadir, subst_model, prefix)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python run_fastrfs.py datadir subst_model")
    sys.exit(1)
  run_fastrfs(sys.argv[1], sys.argv[2])
  



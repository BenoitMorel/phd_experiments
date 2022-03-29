import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/mappings')
sys.path.insert(0, 'tools/families')
import experiments as exp
import time
import saved_metrics
import ete3
import get_dico
import random
import species_analyze
import time
from run_concatenation import build_supermatrix


def exec_superfastme(subst_model, is_dna, cores, run_dir, supermatrix_path, output_tree):
  print("un raxml")  
  command = []
  command.append(exp.fastme_exec)
  command.append("-i")
  command.append(supermatrix_path)
  if (is_dna):
    command.append("-d" + subst_model)
  else:
    command.append("-p" + subst_model)
  command.append("-o")
  command.append(output_tree)
  command.append("--seed")
  command.append("40")
  command.append("--spr")
  #command.append("-T")
  #command.append(str(cores))
  FNULL = open(os.devnull, 'w')
  print("running " + " ".join(command))
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
  stdout, stderr = process.communicate()

def run_superfastme(datadir, concatenation_mode,  subst_model, is_dna, cores, additional_arguments = []):
  run_name = "superfastme-" + concatenation_mode
  run_dir = fam.get_run_dir(datadir, subst_model,  run_name)
  supermatrix_path = os.path.join(run_dir, "supermatrix.fasta")
  partition_path = os.path.join(run_dir, "supermatrix.part")
  shutil.rmtree(run_dir, True)
  os.makedirs(run_dir)
  sites = build_supermatrix(datadir, subst_model, supermatrix_path, partition_path, concatenation_mode, "iphylip_relaxed")
  cores = min(cores, int(sites / 500))
  
  start = time.time()
  
  output_tree = os.path.join(run_dir, "output.newick")
  exec_superfastme(subst_model, is_dna,  cores, run_dir, supermatrix_path, output_tree)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  dest = fam.get_species_tree(datadir, subst_model, run_name)
  shutil.copy(output_tree, dest)

if __name__ == "__main__":
  min_args_number = 6
  if (len(sys.argv) < min_args_number):
    print("syntax: python run_concatenation.py datadir concatenation_mode subst_model is_dna cores")
    print("concatenation modes can be: ")
    print("- min: randomly take ONE gene from each family and each species")
    print("- max: apply min, remove the selected genes, and restart until there is not gene left")
    print("- single: only take single-copy gene families")
    sys.exit(1)
  datadir = sys.argv[1]
  concatenation_mode = sys.argv[2]
  subst_model = sys.argv[3]
  is_dna = sys.argv[4] != "0"
  cores = int(sys.argv[5])
  additional_arguments = sys.argv[min_args_number:]
  assert(concatenation_mode in ["min", "max", "single"])
  run_superfastme(datadir, concatenation_mode, subst_model, cores ,additional_arguments)
  species_analyze.analyze(datadir)


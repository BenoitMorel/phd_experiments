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
import ete3 
import species_analyze
from run_astral_pro import init_gene_trees_file
from run_astral_pro import init_mapping_file


def exec_aster(gene_trees_file, mapping_file, output_species_tree_file, multi, cores):
  command = []
  tmp_output_species_tree_file = output_species_tree_file + ".tmp"
  if (multi):
    command.append(exp.astralpro2_exec)
  else:
    command.append(exp.astral2_exec)
  command.append("-a")
  command.append(mapping_file)
  command.append("-o")
  command.append(tmp_output_species_tree_file)
  command.append("-t")
  command.append(str(cores))
  command.append(gene_trees_file)
  FNULL = open(os.devnull, 'w')
  subprocess.check_call(command)
  shutil.move(tmp_output_species_tree_file, output_species_tree_file)
  

def run_aster(datadir, method, subst_model, multi, cores = 1):
  run_name = "aster_"
  if (multi):
    run_name = "astralpro2_"
  #if (cores > 1):
  #  run_name = run_name + "t" + str(cores) + "_"
  run_name = run_name +  method
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  gene_trees_file = init_gene_trees_file(datadir, method, subst_model, output_dir)
  mapping_file = init_mapping_file(datadir, output_dir)
  print("Start executing aster")
  start = time.time()
  exec_aster(gene_trees_file, mapping_file, fam.get_species_tree(datadir, subst_model, run_name), multi, cores)
  time1 = (time.time() - start)
  print("Runtime: " + str(time1) + "s")
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 

if (__name__ == "__main__"):
  if (len(sys.argv) != 6):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_trees subst_model is_multi cores")
    sys.exit(1)
  datadir = sys.argv[1]
  method = sys.argv[2]
  subst_model = sys.argv[3]
  multi = int(sys.argv[4]) != 0
  cores = int(sys.argv[5])
  run_aster(datadir, method, subst_model, multi, cores)
  species_analyze.analyze(sys.argv[1])
  





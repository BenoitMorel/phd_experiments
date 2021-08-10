import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/raxml/')
import experiments as exp
import raxml_get_tca_score as tca
import time
import saved_metrics
import run_raxml_supportvalues as raxml
import sequence_model

def run_pargenes(datadir, pargenes_dir, subst_model, starting_trees, cores):
  parsimony_trees = int(starting_trees) // 2
  random_trees = starting_trees - parsimony_trees
  raxml_command = ""
  run_modeltest = (subst_model == "bestAA" or subst_model == "bestNT")
  if (not run_modeltest):
    raxml_command +="--model " + sequence_model.get_raxml_model(subst_model) + " --blopt nr_safe"
  command = []
  command.append(exp.python())
  command.append(exp.pargenes_script_debug)
  command.append("-a")
  command.append(os.path.join(datadir, "alignments"))
  command.append("-o")
  command.append(pargenes_dir)
  command.append("-c")
  command.append(str(cores))
  command.append("-s")
  command.append(str(random_trees))
  if (parsimony_trees > 0):
    command.append("-p")
    command.append(str(parsimony_trees))
  if (len(raxml_command) > 0):
    command.append("-R")
    command.append(raxml_command)
  if (run_modeltest):
    command.append("-m")
    if (subst_model == "bestAA"):
      command.append("-d")
      command.append("aa")
  command.append("--continue")
  try:
    subprocess.check_call(command, stdout = sys.stdout)
  except:
    command[0] = exp.python()
    print(" ".join(command))
    subprocess.check_call(command, stdout = sys.stdout)


def export_pargenes_trees(pargenes_dir, starting_trees, subst_model, datadir):
  families_dir = os.path.join(datadir, "families")
  ml_trees_dir = os.path.join(pargenes_dir, "mlsearch_run", "results")
  for family in os.listdir(ml_trees_dir):
    trees_file = os.path.join(ml_trees_dir, family, "sorted_ml_trees.newick")
    family = "_".join(family.split("_")[:-1])
    new_raxml_trees = fam.get_raxml_trees(datadir, starting_trees, subst_model, family)
    try:
      shutil.copyfile(trees_file, new_raxml_trees)
    except:
      print("Cannot copy " + trees_file + " to " + new_raxml_trees)
      pass

def run_pargenes_and_extract_trees(datadir, subst_model, starting_trees, cores, restart = False):
  saved_metrics_key = "raxml-ng" + str(starting_trees)
  pargenes_dir = saved_metrics_key
  pargenes_dir = fam.get_run_dir(datadir, subst_model, pargenes_dir)
  if (not restart):
    shutil.rmtree(pargenes_dir, True)
  start = time.time()
  run_pargenes(datadir, pargenes_dir, subst_model, starting_trees, cores)
  saved_metrics.save_metrics(datadir, fam.get_run_name(saved_metrics_key, subst_model), (time.time() - start), "runtimes") 
  export_pargenes_trees(pargenes_dir, starting_trees, subst_model, datadir)

if __name__ == "__main__":
  if (len(sys.argv) < 6):
    print("syntax: python " + os.path.basename(__file__) + " datadir subst_model starting_trees cores restart")
    sys.exit(1)
  dataset = sys.argv[1]
  subst_model = sys.argv[2]
  starting_trees = int(sys.argv[3])
  cores = int(sys.argv[4])
  restart = int(sys.argv[5]) == 1
  run_pargenes_and_extract_trees(dataset, subst_model, starting_trees, cores, restart = restart)

  


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

def run_pargenes(datadir, pargenes_dir, subst_model, starting_trees, bs_trees, cores):
  
  raxml_command = "--model " + subst_model + " --blopt nr_safe"
  command = []
  command.append("python")
  command.append(exp.pargenes_script_debug)
  command.append("-a")
  command.append(os.path.join(datadir, "alignments"))
  command.append("-b")
  command.append(str(bs_trees))
  command.append("-o")
  command.append(pargenes_dir)
  command.append("-c")
  command.append(str(cores))
  command.append("-s")
  command.append(str(starting_trees))
  command.append("-R")
  command.append(raxml_command)
  command.append("--continue")
  try:
    subprocess.check_call(command, stdout = sys.stdout)
  except:
    command[0] = "python2.7"
    print(" ".join(command))
    subprocess.check_call(command, stdout = sys.stdout)


def export_pargenes_trees(pargenes_dir, subst_model, datadir):
  families_dir = os.path.join(datadir, "families")
  # tca scores
  concatenated_dir = os.path.join(pargenes_dir, "concatenated_bootstraps")
  for concatenation in os.listdir(concatenated_dir):
    family = "_".join(concatenation.split("_")[:-1]) # remove everything after the last
    tca_score = 0.0
    try:
      tca_score = tca.get_tca(os.path.join(concatenated_dir, concatenation), family)
    except:
      print("failed to extract tca score for " + concatenation)
      continue
    try:
      output = os.path.join(families_dir, family, "tca.txt")
      open(output, "w").write(str(tca_score))
    except:
      continue
  # support trees
  support_trees_dir = os.path.join(pargenes_dir, "supports_run", "results")
  for support_tree in os.listdir(support_trees_dir):
    if (not support_tree.endswith("support")):
      continue
    family = "_".join(support_tree.split("_")[:-1]) # remove everything after the last _
    old_raxml_tree = os.path.join(support_trees_dir, support_tree)
    new_raxml_tree = fam.get_raxml_tree(datadir, subst_model, family)
    shutil.copyfile(old_raxml_tree, new_raxml_tree)

  # ml trees
  ml_trees_dir = os.path.join(pargenes_dir, "mlsearch_run", "results")
  for family in os.listdir(ml_trees_dir):
    trees_file = os.path.join(ml_trees_dir, family, "sorted_ml_trees.newick")
    best_model = os.path.join(ml_trees_dir, family, family + ".raxml.bestModel")
    if (not os.path.isfile(trees_file)):
      trees_file = os.path.join(ml_trees_dir, family, family + ".raxml.bestTree")
    family = "_".join(family.split("_")[:-1]) # remove everything after the last _
    new_raxml_tree = fam.get_raxml_multiple_trees(datadir, subst_model, family)
    try:
      shutil.copyfile(trees_file, new_raxml_tree)
    except:
      print("Cannot copy " + trees_file + " to " + new_raxml_tree)
      pass
    new_best_model = fam.get_raxml_best_model(datadir, subst_model, family)
    shutil.copyfile(best_model, new_best_model)
  # clean
  garbage_dir = os.path.join(datadir, "garbage")
  try:
    os.makedirs(garbage_dir)
  except:
    pass
  for family in os.listdir(families_dir):
    if (not os.path.isfile(fam.get_raxml_tree(datadir, subst_model, family))): 
      print("Cleaning family " + family)
      shutil.move(os.path.join(families_dir, family), garbage_dir)

def run_pargenes_and_extract_trees(datadir, subst_model, starting_trees, bs_trees, cores, pargenes_dir = "pargenes", extract_trees = True):
  saved_metrics_key = "RAxML-NG"
  if (pargenes_dir != "pargenes"):
    saved_metrics_key = pargenes_dir
  pargenes_dir = fam.get_run_dir(datadir, subst_model, pargenes_dir)
  shutil.rmtree(pargenes_dir, True)
  start = time.time()
  run_pargenes(datadir, pargenes_dir, subst_model, starting_trees, bs_trees, cores)
  saved_metrics.save_metrics(datadir, saved_metrics_key, (time.time() - start), "runtimes") 
  if (extract_trees):
    export_pargenes_trees(pargenes_dir, subst_model, datadir)

if __name__ == "__main__":
  for dataset in os.listdir(sys.argv[1]):
    dataset = os.path.join(sys.argv[1], dataset)
    print(dataset)
    try:
      run_pargenes_and_extract_trees(dataset, "1", "10", "100", 40)
      print("ok")
    except:
      print("ko")
  exit(0) 
  
  if (len(sys.argv) != 6):
    print("syntax: python run_raxml_supportvalues.py datadir subst_model starting_trees bs_trees cores")
    sys.exit(1)
  dataset = sys.argv[1]
  subst_model = sys.argv[2]
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
  
  run_pargenes_and_extract_trees(dataset, subst_model, starting_trees, bs_trees, cores)

  

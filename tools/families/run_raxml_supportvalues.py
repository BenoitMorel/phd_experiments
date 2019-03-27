import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/raxml/')
import experiments as exp
import raxml_get_tca_score as tca
import time
import runtimes
import run_raxml_supportvalues as raxml

def run_pargenes(dataset_dir, pargenes_dir, is_dna, starting_trees, bs_trees, cores):
  command = []
  command.append("python")
  command.append(exp.pargenes_script)
  command.append("-a")
  command.append(os.path.join(dataset_dir, "alignments"))
  command.append("-b")
  command.append(str(bs_trees))
  command.append("-o")
  command.append(pargenes_dir)
  command.append("-c")
  command.append(str(cores))
  command.append("-s")
  command.append(str(starting_trees))
  if (is_dna):
    command.append("-r")
    command.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "raxml_command.txt"))
  else:
    command.append("-d")
    command.append("aa")
    command.append("-r")
    command.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "raxml_command_prot.txt"))
  command.append("--continue")
  subprocess.check_call(command, stdout = sys.stdout)


def export_pargenes_trees(pargenes_dir, dataset_dir):
  families_dir = os.path.join(dataset_dir, "families")
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
    new_raxml_tree = os.path.join(families_dir, family, "raxmlGeneTree.newick")
    shutil.copyfile(os.path.join(support_trees_dir, support_tree), new_raxml_tree)
  # ml trees
  ml_trees_dir = os.path.join(pargenes_dir, "mlsearch_run", "results")
  for family in os.listdir(ml_trees_dir):
    trees_file = os.path.join(ml_trees_dir, family, "sorted_ml_trees.newick")
    best_model = os.path.join(ml_trees_dir, family, family + ".raxml.bestModel")
    if (not os.path.isfile(trees_file)):
      trees_file = os.path.join(ml_trees_dir, family, family + ".raxml.bestTree")
    family = "_".join(family.split("_")[:-1]) # remove everything after the last _
    new_raxml_tree = os.path.join(families_dir, family, "raxmlGeneTrees.newick")
    try:
      shutil.copyfile(trees_file, new_raxml_tree)
    except:
      print("Cannot copy " + trees_file + " to " + new_raxml_tree)
      pass
    new_best_model = os.path.join(families_dir, family, "raxmlBestModel.txt")
    shutil.copyfile(best_model, new_best_model)
  # clean
  garbage_dir = os.path.join(dataset_dir, "garbage")
  try:
    os.makedirs(garbage_dir)
  except:
    pass
  for family in os.listdir(families_dir):
    if (not os.path.isfile(os.path.join(families_dir, family, "raxmlGeneTree.newick"))): 
      print("Cleaning family " + family)
      shutil.move(os.path.join(families_dir, family), garbage_dir)

def run_pargenes_and_extract_trees(dataset_dir, is_dna, starting_trees, bs_trees, cores, pargenes_dir = "pargenes", extract_trees = True):
  runtimes_key = "RAxML-NG"
  if (pargenes_dir != "pargenes"):
    runtimes_key = pargenes_dir
  pargenes_dir = os.path.join(dataset_dir, pargenes_dir)
  start = time.time()
  run_pargenes(dataset_dir, pargenes_dir, is_dna, starting_trees, bs_trees, cores)
  runtimes.save_elapsed_time(dataset_dir, runtimes_key, (time.time() - start)) 
  if (extract_trees):
    export_pargenes_trees(pargenes_dir, dataset_dir)

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
    print("syntax: python run_raxml_supportvalues.py dataset_dir is_dna starting_trees bs_trees cores")
    sys.exit(1)
  dataset = sys.argv[1]
  is_dna = sys.argv[2]
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
  
  run_pargenes_and_extract_trees(dataset, is_dna, starting_trees, bs_trees, cores)

  

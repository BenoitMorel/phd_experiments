import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp

def run_pargenes(dataset_dir, pargenes_dir, starting_trees, bs_trees, cores):
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
  command.append("-R")
  command.append("--model GTR --blopt nr_safe")
  subprocess.check_call(command)


def export_pargenes_trees(pargenes_dir, dataset_dir):
  families_dir = os.path.join(dataset_dir, "families")
  # support trees
  support_trees_dir = os.path.join(pargenes_dir, "supports_run", "results")
  for support_tree in os.listdir(support_trees_dir):
    if (not support_tree.endswith("support")):
      continue
    family = "_".join(support_tree.split("_")[:-1]) # remove everything after the last _
    new_raxml_tree = os.path.join(families_dir, family, "raxmlGeneTree.newick")
    print(os.path.join(support_trees_dir, support_tree))
    print(new_raxml_tree)
    shutil.copyfile(os.path.join(support_trees_dir, support_tree), new_raxml_tree)
  # ml trees
  ml_trees_dir = os.path.join(pargenes_dir, "mlsearch_run", "results")
  for family in os.listdir(ml_trees_dir):
    trees_file = os.path.join(ml_trees_dir, family, "sorted_ml_trees.newick")
    if (not os.path.isfile(trees_file)):
      trees_file = os.path.join(ml_trees_dir, family, family + ".raxml.bestTree")
    family = "_".join(family.split("_")[:-1]) # remove everything after the last _
    new_raxml_tree = os.path.join(families_dir, family, "raxmlGeneTrees.newick")
    shutil.copyfile(trees_file, new_raxml_tree)
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

def run_pargenes_and_extract_trees(dataset_dir, starting_trees, bs_trees, cores):
  pargenes_dir = os.path.join(dataset_dir, "pargenes")
  run_pargenes(dataset_dir, pargenes_dir, starting_trees, bs_trees, cores)
  export_pargenes_trees(pargenes_dir, dataset_dir)

if __name__ == "__main__":
  if (len(sys.argv) != 5):
    print("syntax: python run_raxml_supportvalues.py dataset_dir starting_trees bs_trees cores")
    sys.exit(1)
  dataset = sys.argv[1]
  starting_trees = sys.argv[2]
  bs_trees = sys.argv[3]
  cores = int(sys.argv[4])
  run_pargenes_and_extract_trees(dataset, starting_trees, bs_trees, cores)

  

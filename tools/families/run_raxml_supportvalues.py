import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp

def run_pargenes(dataset_dir, pargenes_dir, bs_trees, cores):
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
  command.append("-R")
  command.append("--model GTR")
  subprocess.check_call(command)


def export_pargenes_trees(pargenes_dir, dataset_dir):
  support_trees_dir = os.path.join(pargenes_dir, "supports_run", "results")
  families_dir = os.path.join(dataset_dir, "families")
  for support_tree in os.listdir(support_trees_dir):
    if (not support_tree.endswith("support")):
      continue
    family = support_tree.split("_")[0] + "_pruned"
    new_raxml_tree = os.path.join(families_dir, family, "raxmlGeneTree.newick")
    shutil.copyfile(os.path.join(support_trees_dir, support_tree), new_raxml_tree)
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

def run_pargenes_and_extract_trees(dataset_dir, bs_trees, cores):
  pargenes_dir = os.path.join(dataset_dir, "pargenes")
  run_pargenes(dataset_dir, pargenes_dir, bs_trees, cores)
  export_pargenes_trees(pargenes_dir, dataset_dir)

if __name__ == "__main__":
  if (len(sys.argv) != 4):
    print("syntax: python run_raxml_supportvalues.py dataset_dir bs_trees cores")
    sys.exit(1)
  dataset = sys.argv[1]
  bs_trees = sys.argv[2]
  cores = int(sys.argv[3])
  run_pargenes_and_extract_trees(dataset, bs_trees, cores)

  

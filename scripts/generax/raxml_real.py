import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import saved_metrics
import rf_cells
import experiments as exp
import shutil
import time
import fam
import create_random_tree
import random



def run_raxml_normal(datadir, subst_model, family, cores, output_dir):
  exp.reset_dir(output_dir)
  executable = exp.raxml_exec
  starting_tree = os.path.join(output_dir, "starting_random.newick")
  output_tree = os.path.join(output_dir, "raxml_normal_tree.newick")
  alignment = fam.get_alignment(datadir, family)
  print("Alignment file: " + alignment)
  random.seed(1234)
  create_random_tree.create_random_tree(alignment, starting_tree)
  commands = []
  commands.append(executable)
  commands.append("--msa")
  commands.append(alignment)
  commands.append("--prefix")
  commands.append(os.path.join(output_dir, family))
  commands.append("--model")
  commands.append(subst_model)
  commands.append("--tree")
  commands.append(starting_tree)
  commands.append("--seed")
  commands.append("45")
  print("Executing " + " ".join(commands))
  subprocess.check_call(commands)

if (__name__ == "__main__"): 
  
  if (len(sys.argv) != 6):
    print("Syntax: datadir subst_model family cores output_dir")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  family = sys.argv[3]
  cores = sys.argv[4]
  output_dir = sys.argv[5]
  run_raxml_normal(datadir, subst_model, family, cores, output_dir) 



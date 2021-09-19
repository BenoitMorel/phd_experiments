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


def run_raxml_light(datadir, subst_model, family, cores, output_dir):
  exp.reset_dir(output_dir)
  executable = exp.generax_exec
  starting_tree = os.path.join(output_dir, "starting_random.newick")
  output_tree = os.path.join(output_dir, "raxml_light_tree.newick")
  output_model = os.path.join(output_dir, "raxml_light_model.txt")
  output_stats = os.path.join(output_dir, "raxml_light_stats.txt")
  alignment = fam.get_alignment(datadir, family)
  print("Alignment file: " + alignment)
  random.seed(1234);
  create_random_tree.create_random_tree(alignment, starting_tree)
  commands = []
  commands.append(executable)
  commands.append("raxmlLight")
  commands.append(starting_tree)
  commands.append(alignment)
  commands.append(subst_model)
  commands.append(output_tree)
  commands.append(output_model)
  commands.append(output_stats)
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
  run_raxml_light(datadir, subst_model, family, cores, output_dir) 



import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam


def run_genetegrator_bench(dataset, starting_species_tree, gene_trees, subst_model, cores, additional_args = []):
  command = []
  command.append(exp.python())
  command.append(os.path.join(exp.scripts_root, "generax/launch_tegrator.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append(starting_species_tree)
  command.append(gene_trees)
  command.append("normald")
  command.append(str(cores))
  for arg in additional_args:
    command.append(arg)
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)






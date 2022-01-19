import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam


def run_asteroid(dataset, gene_trees, subst_model, cores, additional_args = []):
  command = []
  command.append(exp.python())
  command.append(os.path.join(exp.scripts_root, "asteroid/launch_asteroid.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append(starting_species_tree)
  command.append(gene_trees)
  command.append("normald")
  command.append(str(cores))
  command.append("--rec-model")
  command.append(rec_model)
  if (prune):
    command.append("--prune-species-tree")
  command.append("--per-family-rates")
  command.append("--skip-family-filtering")
  command.append("--do-not-reconcile")
  for arg in additional_args:
    command.append(arg)
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)


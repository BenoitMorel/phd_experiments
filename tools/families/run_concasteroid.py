import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam


def run_concasteroid(dataset, subst_model, is_dna, cores, additional_args = []):
  command = []
  command.append(exp.python())
  command.append(os.path.join(exp.scripts_root, "asteroid/launch_concasteroid.py"))
  command.append(dataset)
  command.append(subst_model)
  if (is_dna):
    command.append("1")
  else:
    command.append("0")
  command.append(str(cores))
  for arg in additional_args:
    command.append(arg)
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)



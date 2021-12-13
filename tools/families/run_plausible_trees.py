import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/raxml/')
import experiments as exp
import time
import saved_metrics
import run_raxml_supportvalues as raxml
import sequence_model


def run_pargenes(datadir, pargenes_dir, subst_model,  starting_trees, cores):
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
  command.append(str(starting_trees / 2))
  command.append("-p")
  command.append(str(starting_trees - starting_trees / 2))
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


def export_pargenes_trees(pargenes_dir, subst_model, samples, datadir):
  families_dir = os.path.join(datadir, "families")
  # tca scores
  concatenated_dir = os.path.join(pargenes_dir, "concatenated_bootstraps")
  if (os.path.isdir(concatenated_dir)):
    for concatenation in os.listdir(concatenated_dir):
      family = "_".join(concatenation.split("_")[:-1]) # remove everything after the last
      src = os.path.join(concatenated_dir, concatenation)
      dest = fam.get_plausible_trees(datadir, samples, subst_model, family)
      shutil.copyfile(src, dest)



def run_pargenes_and_extract_trees(datadir, subst_model, samples, cores, pargenes_dir = "plausible", extract_trees = True, restart = False):


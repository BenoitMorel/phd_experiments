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


def run_pargenes(datadir, pargenes_dir, subst_model,  samples, cores):
  raxml_command = ""
  run_modeltest = (subst_model == "bestAA" or subst_model == "bestNT")
  if (not run_modeltest):
    raxml_command +="--model " + sequence_model.get_raxml_model(subst_model) + " --blopt nr_safe"
  command = []
  command.append(exp.python())
  command.append(exp.pargenes_script_debug)
  command.append("-a")
  command.append(os.path.join(datadir, "alignments"))
  command.append("-b")
  command.append(str(samples))
  command.append("-o")
  command.append(pargenes_dir)
  command.append("-c")
  command.append(str(cores))
  command.append("-s")
  command.append("0")
  command.append("-p")
  command.append("0")
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
      dest = fam.get_bootstrap_trees(datadir, samples, subst_model, family)
      shutil.copyfile(src, dest)



def run_pargenes_and_extract_trees(datadir, subst_model, samples, cores, pargenes_dir = "bootstrap", extract_trees = True, restart = False):
  saved_metrics_key = "bootstrap" + str(samples)
  if (pargenes_dir != "pargenes"):
    saved_metrics_key = pargenes_dir
  print(datadir)
  print(subst_model)
  print(pargenes_dir)
  pargenes_dir = fam.get_run_dir(datadir, subst_model, pargenes_dir)
  if (not restart):
    shutil.rmtree(pargenes_dir, True)
  start = time.time()
  run_pargenes(datadir, pargenes_dir, subst_model,  samples, cores)
  saved_metrics.save_metrics(datadir, fam.get_run_name(saved_metrics_key, subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(os.path.join(pargenes_dir, "mlsearch_run"))
  saved_metrics.save_metrics(datadir, fam.get_run_name(saved_metrics_key, subst_model), (time.time() - start) * lb, "seqtimes") 
  if (extract_trees):
    export_pargenes_trees(pargenes_dir, subst_model, samples, datadir)
  cleanup = True
  if (cleanup):
    shutil.rmtree(pargenes_dir, True)
    
if __name__ == "__main__":
  if (len(sys.argv) < 6):
    print("syntax: python run_raxml_supportvalues.py datadir subst_model samples cores restart")
    sys.exit(1)
  dataset = sys.argv[1]
  subst_model = sys.argv[2]
  samples = int(sys.argv[3])
  cores = int(sys.argv[4])
  restart = int(sys.argv[5]) == 1
  run_pargenes_and_extract_trees(dataset, subst_model, samples, cores, restart = restart)

  


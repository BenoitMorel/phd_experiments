import os
import sys
import subprocess
import shutil
import time
import concurrent.futures
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import saved_metrics
import experiments as exp
import msa_converter
import rename_leaves


def get_mrbayes_output_dir(datadir, subst_model):
  return fam.get_run_dir(datadir, subst_model, "mrbayes_run")

def generate_mrbayes_commands_file(datadir, generations, frequency, runs, chains, subst_model, cores, output_dir):
  output_dir = get_exabayes_output_dir(datadir, subst_model)
  
def extract_mrbayes_results(datadir, burnin):
  output_dir = get_exabayes_output_dir(datadir, subst_model)
  print("TODO IMPLEMENT")

def run_mrbayes_on_families(datadir, generations, frequency, runs, chains, burnin, subst_model, cores, redo):
  output_dir = get_exabayes_output_dir(datadir, subst_model)
  if (not redo):
    shutil.rmtree(output_dir, True)
  exp.mkdir(output_dir)
  parameters = os.path.join(output_dir, "parameters.txt")
  open(parameters, "w").write("Parameters: " + datadir + " " + str(generations) + " " + str(frequency) + " " + str(runs) + " " + str(chains) + " " + str(burnin) + " " + str(int(subst_model)) + " " + str(cores) + " " + str(int(redo)))
  scheduler_commands_file = generate_mrbayes_commands_file(datadir, generations, frequency, runs, chains, subst_model, cores, output_dir)
  start = time.time()
  exp.run_with_scheduler(exp.exabayes_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, "mrbayes", (time.time() - start), "runtimes") 
  print("Finished running mrbayes after " + str(time.time() - start) + "s")
  sys.stdout.flush()
  extract_mrbayes_results(datadir, burnin)


if (__name__== "__main__"):
  if len(sys.argv) != 10:
    print("Syntax error: python run_exabayes.py datadir generations frequency runs chains cores burnin subst_model redo.")
    print(len(sys.argv))
    sys.exit(0)

  datadir = sys.argv[1]
  generations = int(sys.argv[2])
  frequency = int(sys.argv[3])
  runs = int(sys.argv[4])
  chains = int(sys.argv[5])
  burnin = int(sys.argv[6])
  subst_model = sys.argv[7]
  cores = int(sys.argv[8])
  redo = int(sys.argv[9] != 0)
  run_mrbayes_on_families(datadir, generations, frequency, runs, chains, burnin, subst_model, cores, redo)


#


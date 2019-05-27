import os
import sys
import subprocess
import shutil
import time
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "scheduler"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import scheduler
import saved_metrics
import experiments as exp
import msa_converter
import cut_node_names 
import re
import  run_exabayes

# ALE
ALE_SAMPLES = 1

# EXABAYES
EXA_TREES = 2000
EXA_FREQ = 1000
NUM_RUNS = 4
NUM_CHAINS = 4
PER_RUN_BURN_IN = 200
EXA_GEN = EXA_TREES * EXA_FREQ


def get_method_name(with_transfers):
  if (with_transfers):
    return "ale-dtl"
  else:
    return "ale-dl"

def get_observe_run_dir(datadir, with_transfers):
  method_name = get_method_name(with_transfers)
  return os.path.join(datadir, "runs", "ALE", method_name + "_observe_run")

def clean_ALE(datadir):
  try:  
    shutil.rmtree(get_observe_run_dir(datadir, True))
  except:
    pass
  try:
    shutil.rmtree(get_observe_run_dir(datadir, False))
  except:
    pass

def generate_ALE_observe_commands_file(datadir, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      tree_list = os.path.join(fam.get_misc_dir(datadir, family), family + ".treelist")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(tree_list)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def generate_ALE_ml_commands_file(datadir, with_transfers, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  speciesTree = fam.get_phyldog_species_tree(datadir)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      ale_object = os.path.join(fam.get_misc_dir(datadir, family), family + ".treelist.ale")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(speciesTree)
      command.append(ale_object)
      command.append("sample=" + str(ALE_SAMPLES))
      command.append("seed=42")
      if (not with_transfers):
        command.append("tau=0.0")
      command.append("separator=_")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file


def force_move(src, dest_dir):
  shutil.move(src, os.path.join(dest_dir, src))

def extract_trees_from_ale_output(ale_output, output_trees):
  lines = open(ale_output).readlines()
  position = 0
  while True:
    if ("reconciled G-s:" in lines[position]):
      break
    position += 1
  position += 2
  with open(output_trees, "w") as writer:
    for i in range(position, position + ALE_SAMPLES):
      # get rid of weird reconciliation format
      new_lines = re.sub("\.T.[^:]*", "", lines[i])
      new_lines = re.sub("\.[0-9\.]*:", ":", new_lines)
      writer.write(new_lines)

def extract_ALE_results(datadir, ALE_run_dir, with_transfers, families_dir):
  method_name = get_method_name(with_transfers)
  for family in fam.get_families_list(datadir):
    family_misc_dir = fam.get_misc_dir(datadir, family)
    family_trees_dir = fam.get_gene_tree_dir(datadir, family)
    prefix =  "phyldogSpeciesTree.newick_" + family + ".treelist.ale"
    prefixed_output_trees = os.path.join(family_misc_dir, method_name + "_samples_prefixed.newick")
    output_trees = fam.get_ale_tree(datadir, family, method_name)
    extract_trees_from_ale_output(prefix + ".uml_rec", prefixed_output_trees) 
    cut_node_names.remove_prefix_from_trees(prefixed_output_trees, output_trees) 
# clean files
    force_move(prefix + ".uTs", family_misc_dir)
    #force_move(prefix + ".ucons_tree", family_misc_dir)
    force_move(prefix + ".uml_rec", family_misc_dir)

def run_ALE_on_families(datadir, with_transfers, cores):
  method_name = get_method_name(with_transfers)
  observe_output_dir = get_observe_run_dir(datadir, with_transfers)
  ml_output_dir = os.path.join(datadir, "runs",  "ALE", method_name + "_ml_run")
  shutil.rmtree(observe_output_dir, True)
  shutil.rmtree(ml_output_dir, True)
  os.makedirs(observe_output_dir)
  os.makedirs(ml_output_dir)
  commands_observe = generate_ALE_observe_commands_file(datadir, cores, observe_output_dir)
  commands_ml = generate_ALE_ml_commands_file(datadir, with_transfers, cores, ml_output_dir)
  start = time.time()
  scheduler.run_scheduler(commands_observe, exp.ale_observe_exec, cores, observe_output_dir, method_name + "_observe_run.logs")
  scheduler.run_scheduler(commands_ml, exp.ale_ml_exec, cores, ml_output_dir, method_name + "_ml_run.logs")
  saved_metrics.save_metrics(datadir, method_name, (time.time() - start), "runtimes") 
  extract_ALE_results(datadir, ml_output_dir, with_transfers, os.path.join(datadir, "families"))

def run_exabayes_and_ALE(datadir, is_dna, cores, runs = NUM_RUNS, chains = NUM_CHAINS, generations = EXA_GEN, frequency = EXA_FREQ, burnin = PER_RUN_BURN_IN, redo_exabayes = False):
  if (runs < 0):
    runs = NUM_RUNS
  if (chains < 0):
    chains = NUM_CHAINS
  cwd = os.getcwd()
  try:
    run_dir = os.path.join(datadir, "runs")
    parameters = os.path.join(run_dir, "parameters.txt")
    open(parameters, "w").write(datadir + " " + str(int(is_dna)) + " " + str(cores) + " " + str(runs) + " " + str(chains) + " " + str(generations) + " " + str(frequency) + " " + str(burnin) + " " + str(int(redo_exabayes)))
        
    datadir = os.path.abspath(datadir)
    os.chdir(run_dir)
    if (run_exabayes):
      run_exabayes.run_exabayes_on_families(datadir, generations, frequency, runs, chains, burnin, is_dna, cores, redo_exabayes)
    run_ALE_on_families(datadir, True, cores)
    run_ALE_on_families(datadir, False, cores)
    run_exabayes.clean_exabayes(datadir)
    clean_ALE(datadir)
  finally:
    os.chdir(cwd)

def restart_exabayes_and_ALE(datadir, cores):
  run_dir = os.path.join(datadir, "runs")
  parameters_file = os.path.join(run_dir, "parameters.txt")
  p = open(parameters_file).readlines()[0].split()
  print("Restarting exabayes and ALE with the parameters: " + " ".join(p))
  datadir = p[0]
  is_dna = int(p[1])
  runs = int(p[3])
  chains =  int(p[4])
  generations = int(p[5])
  frequency = int(p[6])
  burnin = int(p[7])
  redo_exabayes = True
  run_exabayes_and_ALE(datadir, is_dna, cores, runs, chains, generations, frequency, burnin, redo_exabayes)

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_ALE.py datadir is_dna cores.")
    sys.exit(0)

  datadir = sys.argv[1]
  is_dna = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])
  redo_exabayes_and_ALE(datadir, is_dna, cores)

#


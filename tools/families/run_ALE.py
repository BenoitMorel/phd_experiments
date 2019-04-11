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
ALE_SAMPLES = 10

# EXABAYES
EXA_TREES = 2000
EXA_FREQ = 1000
NUM_RUNS = 4
NUM_CHAINS = 4
PER_RUN_BURN_IN = 100
EXA_GEN = EXA_TREES * EXA_FREQ

def generate_ALE_observe_commands_file(dataset_dir, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      tree_list = os.path.join(family_dir, "misc", family + ".treelist")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(tree_list)
      #command.append("burnin=" + str(BURN_IN))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def generate_ALE_ml_commands_file(dataset_dir, with_transfers, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  speciesTree = os.path.join(dataset_dir, "phyldogSpeciesTree.newick")
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      ale_object = os.path.join(family_dir, "misc", family + ".treelist.ale")
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

def extract_ALE_results(dataset_dir, ALE_run_dir, with_transfers, families_dir):
  method_name = "ALE-DL"
  if (with_transfers):
    method_name = "ALE-DTL"
  for family in os.listdir(families_dir):
    family_misc_dir = fam.getMiscDir(dataset_dir, family)
    family_trees_dir = fam.getTreesDir(dataset_dir, family)
    prefix =  "phyldogSpeciesTree.newick_" + family + ".treelist.ale"
    prefixed_output_trees = os.path.join(family_misc_dir, method_name + "_samples_prefixed.newick")
    output_trees = fam.getALETree(dataset_dir, family, method_name)
    extract_trees_from_ale_output(prefix + ".uml_rec", prefixed_output_trees) 
    cut_node_names.remove_prefix_from_trees(prefixed_output_trees, output_trees) 
# clean files
    force_move(prefix + ".uTs", family_misc_dir)
    #force_move(prefix + ".ucons_tree", family_misc_dir)
    force_move(prefix + ".uml_rec", family_misc_dir)

def run_ALE_on_families(dataset_dir, with_transfers, cores):
  method_name = "ALE-DL"
  if (with_transfers):
    method_name = "ALE-DTL"
  observe_output_dir = os.path.join(dataset_dir, "runs",  "ALE", method_name + "_observe_run")
  ml_output_dir = os.path.join(dataset_dir, "runs",  "ALE", method_name + "_ml_run")
  shutil.rmtree(observe_output_dir, True)
  shutil.rmtree(ml_output_dir, True)
  os.makedirs(observe_output_dir)
  os.makedirs(ml_output_dir)
  commands_observe = generate_ALE_observe_commands_file(dataset_dir, cores, observe_output_dir)
  commands_ml = generate_ALE_ml_commands_file(dataset_dir, with_transfers, cores, ml_output_dir)
  start = time.time()
  scheduler.run_scheduler(commands_observe, exp.ale_observe_exec, cores, observe_output_dir, method_name + "_observe_run.logs")
  scheduler.run_scheduler(commands_ml, exp.ale_ml_exec, cores, ml_output_dir, method_name + "_ml_run.logs")
  saved_metrics.save_metrics(dataset_dir, method_name, (time.time() - start), "runtimes") 
  extract_ALE_results(dataset_dir, ml_output_dir, with_transfers, os.path.join(dataset_dir, "families"))

def run_exabayes_and_ALE(dataset_dir, is_dna, cores):
  fam.init_dataset_dir(dataset_dir)
  run_exabayes.run_exabayes_on_families(dataset_dir, EXA_GEN, EXA_FREQ, NUM_RUNS, NUM_CHAINS, PER_RUN_BURN_IN, is_dna, cores)
  run_ALE_on_families(dataset_dir, True, cores)
  run_ALE_on_families(dataset_dir, False, cores)

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_ALE.py dataset_dir is_dna cores.")
    sys.exit(0)

  dataset_dir = sys.argv[1]
  is_dna = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])
  run_exabayes_and_ALE(dataset_dir, is_dna, cores)

#


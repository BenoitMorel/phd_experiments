import os
import sys
import subprocess
import shutil
import time
import utils
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import saved_metrics
import experiments as exp
import msa_converter
import cut_node_names 
import re

SAMPLES_NUMBER = 1
BURN_IN = 100

def generate_ALE_observe_commands_file(dataset_dir, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      ALE_dir = os.path.join(family_dir, "ALE")
      tree_list = os.path.join(family_dir, "misc", family + ".treelist")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(tree_list)
      command.append("burnin=" + str(BURN_IN))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def generate_ALE_ml_commands_file(dataset_dir, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  speciesTree = os.path.join(dataset_dir, "phyldogSpeciesTree.newick")
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      ALE_dir = os.path.join(family_dir, "ALE")
      ale_object = os.path.join(family_dir, "misc", family + ".treelist.ale")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(speciesTree)
      command.append(ale_object)
      command.append("sample=" + str(SAMPLES_NUMBER))
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
  pattern = re.compile("[0-9]*:")
  with open(output_trees, "w") as writer:
    for i in range(position, position + SAMPLES_NUMBER):
      # get rid of weird reconciliation format
      print(lines[i])
      new_lines = re.sub("\.[0-9]*:", ":", lines[i])
      print(new_lines)
      print("")
      writer.write(new_lines)

def extract_ALE_results(ALE_run_dir, families_dir):
  for family in os.listdir(families_dir):
    family_misc_dir = os.path.join(families_dir, family, "misc")
    try:
      os.makedirs(family_misc_dir)
    except:
      pass
    prefix = family + ".treelist.ale"
    prefixed_output_trees = os.path.join(family_misc_dir, "ale_samples_prefixed.newick")
    output_trees = os.path.join(family_misc_dir, "ale_samples.newick")
    extract_trees_from_ale_output(prefix + ".uml_rec", prefixed_output_trees) 
    cut_node_names.remove_prefix_from_trees(prefixed_output_trees, output_trees) 
# clean files
    force_move(prefix + ".uTs", family_misc_dir)
    force_move(prefix + ".ucons_tree", family_misc_dir)
    force_move(prefix + ".uml_rec", family_misc_dir)

def run_ALE_on_families(dataset_dir, cores):
  observe_output_dir = os.path.join(dataset_dir, "ALE_observe_run")
  ml_output_dir = os.path.join(dataset_dir, "ALE_ml_run")
  shutil.rmtree(observe_output_dir, True)
  shutil.rmtree(ml_output_dir, True)
  os.makedirs(observe_output_dir)
  os.makedirs(ml_output_dir)
  commands_observe = generate_ALE_observe_commands_file(dataset_dir, cores, observe_output_dir)
  commands_ml = generate_ALE_ml_commands_file(dataset_dir, cores, ml_output_dir)
  start = time.time()
  utils.run_scheduler(commands_observe, exp.ale_observe_exec, cores, observe_output_dir, "ALE_observe_run.logs")
  utils.run_scheduler(commands_ml, exp.ale_ml_exec, cores, ml_output_dir, "ALE_ml_run.logs")
  saved_metrics.save_metrics(dataset_dir, "ALE", (time.time() - start), "runtimes") 
  extract_ALE_results(ml_output_dir, os.path.join(dataset_dir, "families"))

if (__name__== "__main__"):
  max_args_number = 3
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_ALE.py dataset_dir cores.")
    sys.exit(0)

  dataset_dir = sys.argv[1]
  cores = int(sys.argv[2])
  run_ALE_on_families(dataset_dir, cores)


#


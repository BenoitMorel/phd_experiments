import os
import sys
import subprocess
import shutil
import time
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))
import saved_metrics
import experiments as exp
import msa_converter
import cut_node_names 
import re
import  run_mrbayes

# ALE
ALE_SAMPLES = 100

# EXABAYES
EXA_TREES = 1000000
EXA_FREQ = 1000
NUM_RUNS = 2
NUM_CHAINS = 4
PER_RUN_BURN_IN = 100
EXA_GEN = EXA_TREES * EXA_FREQ


def get_method_name(with_transfers):
  if (with_transfers):
    return "ale-dtl"
  else:
    return "ale-dl"

def get_observe_run_dir(datadir, subst_model, with_transfers):
  return fam.get_run_dir(datadir, subst_model, get_method_name(with_transfers) + "_observe_run")

def get_ml_run_dir(datadir, subst_model, with_transfers):
  return fam.get_run_dir(datadir, subst_model, get_method_name(with_transfers) + "_ml_run")


def clean_ALE(datadir, subst_model):
  try:  
    shutil.rmtree(get_observe_run_dir(datadir, subst_model, True))
  except:
    pass
  try:
    shutil.rmtree(get_observe_run_dir(datadir, subst_model, False))
  except:
    pass

def generate_ALE_observe_commands_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      tree_list = run_mrbayes.get_treelist(datadir, subst_model, family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(tree_list)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def generate_ALE_ml_commands_file(datadir, subst_model, with_transfers, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  speciesTree = fam.get_phyldog_species_tree(datadir)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      ale_object = run_mrbayes.get_treelist(datadir, subst_model, family) + ".ale"
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

def extract_ALE_results(datadir, subst_model, ALE_run_dir, with_transfers, families_dir):
  method_name = get_method_name(with_transfers)
  for family in fam.get_families_list(datadir):
    family_misc_dir = fam.get_family_misc_dir(datadir, family)
    family_trees_dir = fam.get_gene_tree_dir(datadir, family)
    prefix =  "phyldogSpeciesTree.newick_" + family + ".treelist.ale"
    prefixed_output_trees = os.path.join(family_misc_dir, method_name + "_samples_prefixed.newick")
    output_trees = fam.get_ale_tree(datadir, subst_model, family, method_name)
    extract_trees_from_ale_output(prefix + ".uml_rec", prefixed_output_trees) 
    cut_node_names.remove_prefix_from_trees(prefixed_output_trees, output_trees) 
# clean files
    force_move(prefix + ".uTs", family_misc_dir)
    #force_move(prefix + ".ucons_tree", family_misc_dir)
    force_move(prefix + ".uml_rec", family_misc_dir)

def run_ALE_on_families(datadir, subst_model, with_transfers, cores):
  try:
    cwd = os.getcwd()
    method_name = get_method_name(with_transfers)
    observe_output_dir = get_observe_run_dir(datadir, subst_model, with_transfers)
    ml_output_dir = get_ml_run_dir(datadir, subst_model, with_transfers)
    shutil.rmtree(observe_output_dir, True)
    shutil.rmtree(ml_output_dir, True)
    os.makedirs(observe_output_dir)
    os.makedirs(ml_output_dir)
    commands_observe = generate_ALE_observe_commands_file(datadir, subst_model, cores, observe_output_dir)
    commands_ml = generate_ALE_ml_commands_file(datadir, subst_model, with_transfers, cores, ml_output_dir)
    os.chdir(observe_output_dir)
    start = time.time()
    exp.run_with_scheduler(exp.ale_observe_exec, commands_observe, "onecore", cores, observe_output_dir, method_name + "_ml_run.logs")
    time1 = (time.time() - start)
    os.chdir(ml_output_dir)
    start = time.time()
    exp.run_with_scheduler(exp.ale_ml_exec, commands_ml, "onecore", cores, ml_output_dir, method_name + "_ml_run.logs")
    time2 = (time.time() - start)
    cwd = os.getcwd()
    lb1 = fam.get_lb_from_run(observe_output_dir)
    lb2 = fam.get_lb_from_run(ml_output_dir)
    saved_metrics.save_metrics(datadir, fam.get_run_name(method_name, subst_model), time1 + time2, "runtimes") 
    saved_metrics.save_metrics(datadir, fam.get_run_name(method_name, subst_model), time1 * lb1 + time2 * lb2, "seqtimes") 
    extract_ALE_results(datadir, subst_model, ml_output_dir, with_transfers, os.path.join(datadir, "families"))
  finally:
    cwd = os.getcwd()

def run_ALE(datadir, subst_model, cores):
  cwd = os.getcwd()
  try:
    run_dir = os.path.join(datadir, "runs")
    parameters = os.path.join(run_dir, "parameters.txt")
    datadir = os.path.abspath(datadir)
    os.chdir(run_dir)
    run_ALE_on_families(datadir, subst_model, True, cores)
    run_ALE_on_families(datadir, subst_model, False, cores)
    clean_ALE(datadir, subst_model)
  finally:
    os.chdir(cwd)


if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_ALE.py datadir subst_model cores.")
    sys.exit(0)

  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])
  run_ALE(datadir, subst_model, cores)

#


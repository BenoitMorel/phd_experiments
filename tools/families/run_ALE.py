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


def get_run_name(gene_trees, with_transfers, dated):
  res = ""
  if (dated):
    res += "dated_"
  if (with_transfers):
    res +=  "ale-dtl"
  else:
    res += "ale-dl"
  res += "-"
  res += gene_trees
  return res


def clean_ALE(datadir, subst_model, dated):
  try:  
    shutil.rmtree(get_observe_run_dir(datadir, subst_model, True, dated))
  except:
    pass
  try:
    shutil.rmtree(get_observe_run_dir(datadir, subst_model, False, dated))
  except:
    pass

def generate_ALE_observe_commands_file(datadir, gene_trees, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      tree_list = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
      ale_tree_list = os.path.join(results_dir, family + ".newick")
      exp.relative_symlink(tree_list, ale_tree_list)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(ale_tree_list)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def generate_ALE_ml_commands_file(datadir, subst_model, with_transfers, cores, observe_output_dir, output_dir, additional_arguments):
  obs_results_dir = os.path.join(observe_output_dir, "results")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  os.makedirs(results_dir)
  speciesTree = fam.get_phyldog_species_tree(datadir)
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      ale_tree_list = os.path.join(obs_results_dir, family + ".newick")
      ale_object = ale_tree_list + ".ale"
      
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
      command.extend(additional_arguments)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file


def force_move(src, dest_dir):
  shutil.move(src, os.path.join(dest_dir, src))

def extract_trees_from_ale_output(ale_output, output_trees):
  try:
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
  except:
    print("WARNING: can't extract ale tree " + ale_output + " to " + output_trees)


def extract_ALE_results(datadir, gene_trees, subst_model, ALE_run_dir, with_transfers, families_dir, dated):
  run_name = get_run_name(gene_trees, with_transfers, dated)
  print(ALE_run_dir)
  mlstring = "uml"
  if (dated):
    mlstring = "ml"
  for family in fam.get_families_list(datadir):
    family_misc_dir = fam.get_family_misc_dir(datadir, family)
    family_trees_dir = fam.get_gene_tree_dir(datadir, family)
    prefix =  "phyldogSpeciesTree.newick_" + family + ".newick.ale"
    prefix = os.path.join(ALE_run_dir, prefix)
    prefixed_output_trees = os.path.join(family_misc_dir, run_name + "_samples_prefixed.newick")
    output_trees = fam.get_ale_tree(datadir, subst_model, family, run_name)
    ale_output = prefix +"." + mlstring + "_rec"
    are_prefixed = False
    if (are_prefixed):
      extract_trees_from_ale_output(ale_output, prefixed_output_trees) 
    
      cut_node_names.remove_prefix_from_trees(prefixed_output_trees, output_trees) 
      if (dated):
        cut_node_names.cut_keep_first_elems(output_trees, output_trees, "@", 1) 
    else:
      extract_trees_from_ale_output(ale_output, output_trees) 
      if (dated):
        cut_node_names.cut_keep_first_elems(output_trees, output_trees, "@", 1) 

# clean files
    #if (dated):
    #  force_move(prefix + ".Ts", family_misc_dir)
    #else:
    #  force_move(prefix + ".uTs", family_misc_dir)
    #force_move(prefix + ".ucons_tree", family_misc_dir)
    #force_move(prefix + "." + mlstring + "_rec", family_misc_dir)

def run_ALE_on_families(datadir, gene_trees, subst_model, with_transfers, cores, dated = False, additional_arguments = []):
  cwd = os.getcwd()
  #try:
  run_name = get_run_name(gene_trees, with_transfers, dated)
  output_dir = fam.get_run_dir(datadir, subst_model, run_name)
  observe_output_dir = os.path.join(output_dir, "observe")
  ml_output_dir = os.path.join(output_dir, "ml")
  do_run = True
  if (do_run):
    exp.reset_dir(observe_output_dir)
    exp.reset_dir(ml_output_dir)
    commands_observe = generate_ALE_observe_commands_file(datadir, gene_trees, subst_model, cores, observe_output_dir)
    commands_ml = generate_ALE_ml_commands_file(datadir, subst_model, with_transfers, cores, observe_output_dir,  ml_output_dir, additional_arguments)
    #os.chdir(observe_output_dir)
    start = time.time()
    exp.run_with_scheduler(exp.ale_observe_exec, commands_observe, "onecore", cores, observe_output_dir, run_name + "_ml_run.logs")
    time1 = (time.time() - start)
    os.chdir(ml_output_dir)
    start = time.time()
    if (not dated):
      exp.run_with_scheduler(exp.ale_ml_exec, commands_ml, "onecore", cores, ml_output_dir, run_name + "_ml_run.logs")
    else:
      exp.run_with_scheduler(exp.ale_ml_dated_exec, commands_ml, "onecore", cores, ml_output_dir, run_name + "_ml_run.logs")
    time2 = (time.time() - start)
    os.chdir(cwd)
    lb1 = fam.get_lb_from_run(observe_output_dir)
    lb2 = fam.get_lb_from_run(ml_output_dir)
    saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1 + time2, "runtimes") 
    saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1 * lb1 + time2 * lb2, "seqtimes") 
  
  extract_ALE_results(datadir, gene_trees, subst_model, ml_output_dir, with_transfers, os.path.join(datadir, "families"), dated)
  #finally:
  #  cwd = os.getcwd()

def run_ALE(datadir, gene_trees, subst_model, cores, dated = False, additional_arguments = []):
  cwd = os.getcwd()
  try:
    datadir = os.path.abspath(datadir)
    run_ALE_on_families(datadir, gene_trees, subst_model, True, cores, dated, additional_arguments)
    clean_ALE(datadir, subst_model, dated)
  finally:
    os.chdir(cwd)


if (__name__== "__main__"):
  min_args_number = 5
  if len(sys.argv) < min_args_number:
    print("Syntax error: python run_ALE.py datadir gene_trees subst_model cores [dated] [additional_arguments.")
    sys.exit(0)

  datadir = sys.argv[1]
  gene_trees = sys.argv[2]
  subst_model = sys.argv[3]
  cores = int(sys.argv[4])
  dated = False
  additional_arguments = []
  if (len(sys.argv) > 5):
    dated = (int(sys.argv[5]) != 0)
  if (len(sys.argv) > 6):
    additional_arguments = sys.argv[6:]
  run_ALE(datadir, gene_trees, subst_model, cores, dated, additional_arguments)

#


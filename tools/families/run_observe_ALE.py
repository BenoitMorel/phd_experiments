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
import run_ALE
  
 
def extract_observe_ALE(datadir, gene_trees, subst_model, observe_output_dir):
  results_dir = os.path.join(observe_output_dir, "results")
  for family in fam.get_families_list(datadir):
    src = os.path.join(results_dir, family + ".newick.ale")
    dest = fam.build_gene_tree_path(datadir, subst_model, family, "ccp-" + gene_trees)
    shutil.copyfile(src, dest)

def run_observe_ALE(datadir, gene_trees, subst_model, cores):
  run_name = "ale_observe_" + gene_trees
  output_dir = fam.get_run_dir(datadir, subst_model, run_name)
  observe_output_dir = os.path.join(output_dir, "observe")
  exp.reset_dir(observe_output_dir)
  commands = run_ALE.generate_ALE_observe_commands_file(datadir, gene_trees, subst_model, cores, observe_output_dir)
  exp.run_with_scheduler(exp.ale_observe_exec, commands, "onecore", cores, observe_output_dir, run_name + "_ml_run.logs")
  extract_observe_ALE(datadir, gene_trees, subst_model, observe_output_dir)


if (__name__== "__main__"):
  min_args_number = 5
  if len(sys.argv) < min_args_number:
    print("Syntax error: python run_ALE.py datadir gene_trees subst_model cores")
    sys.exit(0)

  datadir = sys.argv[1]
  gene_trees = sys.argv[2]
  subst_model = sys.argv[3]
  cores = int(sys.argv[4])
  run_observe_ALE(datadir, gene_trees, subst_model, cores)

#



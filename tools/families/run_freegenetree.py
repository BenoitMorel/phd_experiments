import os
import sys
import subprocess
import shutil
import time
import saved_metrics
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import fam
import run_raxml_supportvalues as run_pargenes
import sequence_model
import ete3

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      freegenetree_dir = fam.get_family_misc_dir(datadir, family)
      try:
        os.mkdir(freegenetree_dir)
      except:
        pass
      freegenetree_output = os.path.join(freegenetree_dir, "freegenetree_output." + subst_model + ".newick")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      if ("LG" in subst_model):
        command.append("--aa")
      command.append("--output")
      command.append(freegenetree_output)
      command.append("--msa")
      command.append(fam.get_alignment(datadir, family))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_freegenetree_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    freegenetreeTree = os.path.join(families_dir, family, "misc", "freegenetree_output." + subst_model + ".newick")
    if (os.path.isfile(freegenetreeTree)):
      tree = fam.build_gene_tree_path(datadir, subst_model, family, "freegenetree")
      shutil.copyfile(freegenetreeTree, tree)

def run_freegenetree_on_families(datadir, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "freegenetree_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir)
  
  start = time.time()
  exp.run_with_scheduler(exp.freegenetree_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("freegenetree", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("freegenetree", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_freegenetree_trees(datadir, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_freegenetree datadir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])

  run_freegenetree_on_families(datadir, subst_model, cores)




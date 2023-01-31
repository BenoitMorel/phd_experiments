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
import sequence_model
import ete3

def generate_scheduler_commands_file(datadir, subst_model, samples, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  exp.reset_dir(results_dir)
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("-s")
      command.append(fam.get_alignment(datadir, family))
      command.append("-m")
      command.append(subst_model)
      command.append("-B")
      command.append(str(samples))
      command.append("-pre")
      command.append(os.path.join(results_dir, family))
      command.append("--boot-trees")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
    
def extract_trees(datadir, output_dir, subst_model, samples):
  results_dir = os.path.join(output_dir, "results")
  for family in fam.get_families_list(datadir):
    src = os.path.join(results_dir, family + ".ufboot")
    dest = fam.build_gene_tree_path(datadir, subst_model, family, "ufboot" + str(samples))
    shutil.copyfile(src, dest)

def run_ufboot_on_families(datadir, subst_model, samples, cores):
  run_name = "ufboot" + str(samples)
  output_dir = fam.get_run_dir(datadir, subst_model, run_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, samples, cores, output_dir)
  start = time.time()
  exp.run_with_scheduler(exp.iqtree_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("run_nae", subst_model), (time.time() - start), "runtimes") 
  extract_trees(datadir, output_dir, subst_model, samples)

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir subst_model samples cores.")
    sys.exit(0)

  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  samples = int(sys.argv[3])
  cores = int(sys.argv[4])
  run_ufboot_on_families(datadir, subst_model, samples, cores)





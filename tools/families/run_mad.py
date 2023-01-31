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




def generate_scheduler_commands_file(datadir, gene_trees, subst_model, cores, outputdir):
  results_dir = os.path.join(outputdir, "results")
  scheduler_commands_file = os.path.join(outputdir, "commands.txt")
  
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      gene_tree = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees) 
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(gene_tree)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_mad_trees(datadir, outputdir, gene_trees, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    src = os.path.join(outputdir, "per_job_logs", family + "_out.txt")
    dest = fam.build_gene_tree_path(datadir, subst_model, family, "mad-" + gene_trees) 
    shutil.copy(src, dest)

def run_mad_on_families(datadir, gene_trees, subst_model, cores):
  key = "mad-" + gene_trees
  outputdir = fam.get_run_dir(datadir, subst_model, key + "_run")
  shutil.rmtree(outputdir, True)
  os.makedirs(outputdir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, gene_trees, subst_model, cores, outputdir)
  start = time.time()
  exp.run_with_scheduler(exp.mad_exec, scheduler_commands_file, "onecore", cores, outputdir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name(key, subst_model), (time.time() - start), "runtimes") 
  extract_mad_trees(datadir, outputdir, gene_trees, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir gene_trees subst_model cores.")
    sys.exit(0)


  datadir = sys.argv[1]
  gene_trees = sys.argv[2]
  subst_model = sys.argv[3]
  cores = int(sys.argv[4])

  run_mad_on_families(datadir, gene_trees, subst_model, cores)




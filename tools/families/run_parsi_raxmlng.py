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

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      parsimony_dir = fam.get_family_misc_dir(datadir, family)
      exp.mkdir(parsimony_dir)
      parsimony_output = os.path.join(parsimony_dir, "parsimony_output." + subst_model)
      alignment = os.path.abspath(os.path.join(family_dir, "alignment.msa"))
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(exp.raxml_nompi_exec)
      command.append("--seed")
      command.append("42")
      command.append("--model")
      command.append(subst_model)
      command.append("--msa")
      command.append(alignment)
      command.append("--prefix")
      command.append(parsimony_output)
      command.append("--start")
      command.append("--tree")
      command.append("pars{10}")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def run_parsi_on_families(datadir, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "parsimony_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir)
  
  start = time.time()
  exp.run_with_scheduler(exp.raxml_nompi_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("parsimony", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("parsimony", subst_model), (time.time() - start) * lb, "seqtimes") 
  #extract_parsimony_trees(datadir, subst_model)

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_parsimony datadir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])
  run_parsi_on_families(datadir, subst_model, cores)




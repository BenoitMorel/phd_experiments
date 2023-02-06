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

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir, likelihoods):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      plausiblerax_dir = fam.get_family_misc_dir(datadir, family)
      exp.mkdir(plausiblerax_dir)
      exp.mkdir(fam.get_likelihoods_dir(datadir, family))
      plausiblerax_output = os.path.join(plausiblerax_dir, "plausiblerax_output." + subst_model + ".newick")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      #command.append(exp.plausiblerax_exec)
      #../FastTree/FastTree -out out.newick -nt -gtr
      raxml_tree = fam.get_raxml_tree(datadir, subst_model, family)
      command.append(raxml_tree)
      command.append(fam.get_alignment(datadir, family))
      command.append(subst_model)
      command.append(fam.build_gene_tree_path(datadir, subst_model, family, "plausiblerax"))
      command.append(fam.build_likelihood_path(datadir, subst_model, family, "plausiblerax"))
      if (likelihoods):
        command.append("--likelihoods")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def run_plausiblerax_on_families(datadir, subst_model, cores, likelihoods = False):
  output_dir = fam.get_run_dir(datadir, subst_model, "plausiblerax_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir, likelihoods)
  start = time.time()
  exp.run_with_scheduler(exp.plausiblerax_exec, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("plausiblerax", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("plausiblerax", subst_model), (time.time() - start) * lb, "seqtimes") 
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir subst_model cores.")
    sys.exit(0)

  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])
  run_plausiblerax_on_families(datadir, subst_model, cores)




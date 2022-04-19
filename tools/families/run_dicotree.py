import os
import sys
import shutil
import time
import saved_metrics
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import fam
import sequence_model
import ete3

def generate_scheduler_commands_file(datadir, subst_model, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      dicotree_dir = os.path.join(output_dir, "dico_" + family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(exp.dicotree_script)
      command.append(fam.get_alignment(datadir, family))
      command.append(subst_model)
      command.append(dicotree_dir)
      command.append("20")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_dicotree_trees(datadir, output_dir, subst_model):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    dicotree_dir = os.path.join(output_dir, "dico_" + family)
    dicoasteroid = fam.build_gene_tree_path(datadir, subst_model, family, "dicotree-asteroid-nobl")
    dicoasteroidbl = fam.build_gene_tree_path(datadir, subst_model, family, "dicotree-asteroid-bl")
    dicoastral = fam.build_gene_tree_path(datadir, subst_model, family, "dicotree-astral")
    dicomrptnt = fam.build_gene_tree_path(datadir, subst_model, family, "dicotree-mrp-tnt")
    shutil.copyfile(os.path.join(dicotree_dir, "final_asteroid.newick"), dicoasteroid)
    shutil.copyfile(os.path.join(dicotree_dir, "final_asteroid-bl.newick"), dicoasteroidbl)
    shutil.copyfile(os.path.join(dicotree_dir, "final_mrp_tnt.newick"), dicomrptnt)
    #shutil.copyfile(os.path.join(dicotree_dir, "final_astral.newick"), dicoastral)


def run_dicotree_on_families(datadir, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, subst_model, "dicotree_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, cores, output_dir)
  
  start = time.time()
  exp.run_with_scheduler(exp.python(), scheduler_commands_file, "fork", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, fam.get_run_name("dicotree", subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  saved_metrics.save_metrics(datadir, fam.get_run_name("dicotree", subst_model), (time.time() - start) * lb, "seqtimes") 
  extract_dicotree_trees(datadir, output_dir, subst_model)
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])

  run_dicotree_on_families(datadir, subst_model, cores)




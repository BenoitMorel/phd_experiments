import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam


def run_speciesrax_bench(dataset, rec_model, prune, starting_species_tree, gene_trees, subst_model, cores, additional_args = []):
  command = []
  command.append(exp.python())
  command.append(os.path.join(exp.scripts_root, "generax/launch_speciesrax.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append(starting_species_tree)
  command.append(gene_trees)
  command.append("normald")
  command.append(str(cores))
  command.append("--rec-model")
  command.append(rec_model)
  if (prune):
    command.append("--prune-species-tree")
  command.append("--per-family-rates")
  command.append("--skip-family-filtering")
  command.append("--do-not-reconcile")
  for arg in additional_args:
    command.append(arg)
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)


def run_speciesrax_instance(dataset, starting_tree, with_transfers, run_name, subst_model, strategy, cores = 40, additional_args = []):
  command = []
  command.append(exp.python())
  command.append(os.path.join(exp.scripts_root, "generax/launch_speciesrax.py"))
  command.append(dataset)
  command.append(starting_tree)
  command.append(subst_model)
  command.append("MiniNJ")
  #command.append("random")
  command.append("normald")
  command.append(str(cores))
  if (with_transfers == 2):
    command.append("--rec-model")
    command.append("ParsimonyDL")
  elif (with_transfers == True):
    command.append("--rec-model")
    command.append("UndatedDTL")
  else:
    command.append("--rec-model")
    command.append("UndatedDL")
  command.append("--si-strategy")
  command.append(strategy)
  for arg in additional_args:
    command.append(arg)
  command.append("--skip-family-filtering")
  command.append("--run")
  command.append(fam.get_run_name(run_name, subst_model))
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)
    
  
def run_speciesrax_on_families(dataset_dir, gene_trees, subst_model, cores, transfers = True, strategy = "HYBRID", rates_per_family = True):
  dataset = os.path.basename(dataset_dir)
  run_name = "speciesrax-"
  args = []
  if (transfers):
    run_name += "dtl-"
  else:
    run_name += "dl-"
  run_name += gene_trees + "-"
  if (rates_per_family):
    args.append("--per-family-rates")
    run_name += "perfam-"
  run_name += strategy
  run_speciesrax_instance(dataset, gene_trees, transfers, run_name, subst_model, strategy, cores = 40, additional_args = args)


if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_speciesrax.py dataset_dir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  subst_model = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])
  run_speciesrax_on_families(dataset_dir, subst_model, cores)





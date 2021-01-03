import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam

def run_generax_instance(dataset, starting_tree, with_transfers, method, subst_model, per_sp_rates, optimize_species, cores = 40):
  command = []
  command.append("python")
  if (optimize_species):
    command.append(os.path.join(exp.scripts_root, "generax/launch_speciesrax.py"))
  else:
    command.append(os.path.join(exp.scripts_root, "generax/launch_generax.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append("SPR")
  if (optimize_species):
    command.append("random")
  else:
    command.append("true")
  command.append(starting_tree)
  command.append("normal")
  command.append(str(cores))
  if (with_transfers):
    command.append("--rec-model")
    command.append("UndatedDTL")
  else:
    command.append("--rec-model")
    command.append("UndatedDL")
  if (optimize_species):
    command.append("--optimize-species-tree")
  if (per_sp_rates):
    command.append("--per-species-rates")
    method = method + "-psr"
  command.append("--analyze")
  command.append("no")
  command.append("--run")
  command.append(fam.get_run_name(method, subst_model))
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)
    
  
def run_generax_on_families(dataset_dir, subst_model, cores, raxml = True, random = True, dl = True, dtl = True, optimize_species = False):
  dataset = os.path.basename(dataset_dir)
  key = "generax"
  if (optimize_species):
    key = "speciesrax"
  if (raxml):
    if (dl):
      run_generax_instance(dataset, "raxml-ng", False, key + "-dl-raxml", subst_model, False, optimize_species, cores)
    if (dtl):
      run_generax_instance(dataset, "raxml-ng", True, key + "-dtl-raxml", subst_model, False, optimize_species, cores)
  if (random):
    if (dl):
      run_generax_instance(dataset, "random", False, key + "-dl-random", subst_model, False, optimize_species, cores)
    if (dtl):
      run_generax_instance(dataset, "random", True, key + "-dtl-random", subst_model, False, optimize_species, cores)
  

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_generax.py dataset_dir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  subst_model = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])
  run_generax_on_families(dataset_dir, subst_model, cores)



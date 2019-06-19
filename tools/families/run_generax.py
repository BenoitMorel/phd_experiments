import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam

def run_generax_instance(dataset, starting_tree, with_transfers, method, subst_model, per_sp_rates, cores = 40):
  command = []
  command.append("python")
  command.append(os.path.join(exp.scripts_root, "generax/launch_generax.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append("SPR")
  command.append(starting_tree)
  command.append("normal")
  command.append(str(cores))
  if (with_transfers):
    command.append("--rec-model")
    command.append("UndatedDTL")
  else:
    command.append("--rec-model")
    command.append("UndatedDL")

  if (per_sp_rates):
    command.append("--per-species-rates")
    method = method + "-psr"
  command.append("--run")
  command.append(fam.get_run_name(method, subst_model))
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)
    
  
def run_generax_on_families(dataset_dir, subst_model, cores, raxml = True, random = True, dl = True, dtl = True):
  dataset = os.path.basename(dataset_dir)
  if (raxml):
    if (dl):
      run_generax_instance(dataset, fam.get_run_name("raxml-ng", subst_model), False, "generax-dl-raxml", subst_model, False, cores)
    if (dtl):
      run_generax_instance(dataset, fam.get_run_name("raxml-ng", subst_model), True, "generax-dtl-raxml", subst_model, False, cores)
  if (random):
    if (dl):
      run_generax_instance(dataset, "random", False, "generax-dl-random", subst_model, False, cores)
    if (dtl):
      run_generax_instance(dataset, "random", True, "generax-dtl-random", subst_model, False, cores)
  #run_generax_instance(dataset, "raxml-ng", False, "generax-dl-raxml", subst_model, True, cores)
  #run_generax_instance(dataset, "raxml-ng", True, "generax-dtl-raxml", subst_model, True, cores)
  

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



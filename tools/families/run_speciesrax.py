import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam

def run_speciesrax_instance(dataset, starting_tree, with_transfers, run_name, subst_model, slow, strategy, cores = 40, additional_args = []):
  command = []
  command.append("python")
  command.append(os.path.join(exp.scripts_root, "generax/launch_speciesrax.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append("NJst")
  command.append(starting_tree)
  command.append("normal")
  command.append(str(cores))
  if (with_transfers):
    command.append("--rec-model")
    command.append("UndatedDTL")
  else:
    command.append("--rec-model")
    command.append("UndatedDL")
  if (slow):
    command.append("--species-slow-radius")
    command.append("1")
  command.append("--species-strategy")
  command.append(strategy)
  for arg in additional_args:
    command.append(arg)
  run_name += "-" + strategy
  command.append("--run")
  command.append(fam.get_run_name(run_name, subst_model))
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)
    
  
def run_speciesrax_on_families(dataset_dir, subst_model, cores, dl = True, dtl = True, slow = False, strategy = "SPR"):
  dataset = os.path.basename(dataset_dir)
  if (slow):
    if (dl):
      run_speciesrax_instance(dataset, "raxml-ng", False, "speciesrax-dl-raxml-slow", subst_model, True, strategy, cores)
    if (dtl):
      run_speciesrax_instance(dataset, "raxml-ng", True, "speciesrax-dtl-raxml-slow", subst_model, True, strategy, cores)
  else:
    if (dl):
      run_speciesrax_instance(dataset, "raxml-ng", False, "speciesrax-dl-raxml", subst_model, False, strategy, cores)
    if (dtl):
      run_speciesrax_instance(dataset, "raxml-ng", True, "speciesrax-dtl-raxml", subst_model, False, strategy, cores)


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




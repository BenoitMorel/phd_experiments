import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam

def run_speciesrax_instance(dataset, starting_tree, with_transfers, method, subst_model, cores = 40):
  command = []
  command.append("python")
  command.append(os.path.join(exp.scripts_root, "speciesrax/launch_speciesrax.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append("random")
  command.append(starting_tree)
  command.append("normal")
  command.append(str(cores))
  if (with_transfers):
    command.append("--rec-model")
    command.append("UndatedDTL")
  else:
    command.append("--rec-model")
    command.append("UndatedDL")
  command.append("--run")
  command.append(fam.get_run_name(method, subst_model))
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)
    
  
def run_speciesrax_on_families(dataset_dir, subst_model, cores, dl = True, dtl = True):
  dataset = os.path.basename(dataset_dir)
  if (dl):
    run_speciesrax_instance(dataset, "raxml-ng", False, "speciesrax-dl-raxml", subst_model, cores)
  if (dtl):
    run_speciesrax_instance(dataset, "raxml-ng", True, "speciesrax-dtl-raxml", subst_model, cores)


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




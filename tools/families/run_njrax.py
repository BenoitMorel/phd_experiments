import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import fam

def run_njrax(dataset, nj_strategy, subst_model, cores = 40):
  command = []
  command.append("python")
  command.append(os.path.join(exp.scripts_root, "generax/launch_speciesrax.py"))
  command.append(os.path.basename(os.path.normpath(dataset)))
  command.append(subst_model)
  command.append(nj_strategy)
  command.append("raxml-ng")
  command.append("normal")
  command.append(str(cores))
  command.append("--species-fast-radius")
  command.append("0")
  command.append("--species-strategy")
  command.append("SPR")
  command.append("--run")
  run_name = "njrax-" + nj_strategy
  command.append(fam.get_run_name(run_name, subst_model))
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)
    
  


if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_speciesrax.py dataset_dir subst_model cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset = sys.argv[1]
  nj_strategy = sys.argv[2]
  subst_model = sys.argv[3] 
  cores = int(sys.argv[4])
  run_njrax(dataset, nj_strategy, subst_model, cores)





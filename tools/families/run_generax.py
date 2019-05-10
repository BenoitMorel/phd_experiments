import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp

def run_generax_instance(dataset, starting_tree, with_transfers, run_name, is_dna, cores = 40):
  command = []
  command.append("python")
  command.append(os.path.join(exp.scripts_root, "generax/launch_generax.py"))
  command.append(dataset)
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

  command.append("--run")
  command.append(run_name)
  if (not is_dna):
    command.append("--protein")
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)
    
  
def run_generax_on_families(dataset_dir, is_dna, cores, raxml = True, random = True, dl = True, dtl = True):
  dataset = os.path.basename(dataset_dir)
  if (raxml):
    if (dl):
      run_generax_instance(dataset, "raxml-ng", False, "generax-dl-raxml", is_dna, cores)
    if (dtl):
      run_generax_instance(dataset, "raxml-ng", True, "generax-dtl-raxml", is_dna, cores)
  if (random):
    if (dl):
      run_generax_instance(dataset, "random", False, "generax-dl-random", is_dna, cores)
    if (dtl):
      run_generax_instance(dataset, "random", True, "generax-dtl-random", is_dna, cores)

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_generax.py dataset_dir is_dna cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  is_dna = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])
  run_generax_on_families(dataset_dir, is_dna, cores)



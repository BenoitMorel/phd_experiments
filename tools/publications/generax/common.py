import subprocess
import sys
import os 

sys.path.insert(0, 'scripts')
import experiments as exp

sys.path.insert(0, os.path.join("tools", "families"))
import generate_families_with_jprime as jprime
import run_all

# couples of (species interval, seed) to get a given number of species node
jsim_species_to_params = {}
#jsim_species_to_params[8] = (9, 2)
#jsim_species_to_params[12] = (3, 9)
#jsim_species_to_params[13] = (8, 42)
#jsim_species_to_params[16] = (3, 6)
#jsim_species_to_params[21] = (7, 41)
#jsim_species_to_params[34] = (8, 30)

jsim_species_to_params[5] = (3, 2)
jsim_species_to_params[10] = (3, 13)
jsim_species_to_params[12] = (3, 9)
jsim_species_to_params[14] = (3, 17)
jsim_species_to_params[19] = (3, 11)
jsim_species_to_params[27] = (3, 12)
jsim_species_to_params[41] = (3, 4)

jsim_species_to_params[16] = (3, 6)



protein_datasets = ["swiss", "cyano_simulated", "sub_t0.05_s0.5_cyano_empirical", "sub_t0.01_s0.2_ensembl_8880_15"]

def get_param_from_dataset_name(parameter, dataset):
  if (parameter == "species"):
    return dataset.split("_")[1][1:]
  if (parameter == "families"):
    return dataset.split("_")[2][1:]
  elif (parameter == "sites"):
    return dataset.split("_")[3][5:]
  elif (parameter == "model"):
    return dataset.split("_")[4]
  elif (parameter == "bl"):
    return dataset.split("_")[5][2:]
  elif (parameter == "dup_rate"):
    return dataset.split("_")[6][1:]
  elif (parameter == "loss_rate"):
    return dataset.split("_")[7][1:]
  elif (parameter == "transfer_rate"):
    return dataset.split("_")[8][1:]
  elif (parameter == "tl_ratio"):
    t = get_param_from_dataset_name("transfer_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    return  str(float(t)/float(l))
  elif (parameter == "dl_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    return str(float(d)/float(l))
    return res
  else:
    return "invalid"


def run_generax(dataset, starting_tree, with_transfers, run_name, is_dna):
  command = []
  command.append("python")
  command.append("/home/morelbt/github/phd_experiments/scripts/jointsearch/launch_generax.py")
  command.append(dataset)
  command.append("SPR")
  command.append(starting_tree)
  command.append("normal")
  command.append("40")
  if (with_transfers):
    command.append("--rec-model")
    command.append("UndatedDTL")
  command.append("--run")
  command.append(run_name)
  if (not is_dna):
    command.append("--protein")
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)


def generate_dataset(dataset):
  species = get_param_from_dataset_name("species", dataset)
  species_internal, seed = jsim_species_to_params[int(species)]
  families = int(get_param_from_dataset_name("families", dataset))
  sites = get_param_from_dataset_name("sites", dataset)
  model = get_param_from_dataset_name("model", dataset)
  bl_factor = float(get_param_from_dataset_name("bl", dataset))
  d = get_param_from_dataset_name("dup_rate", dataset)
  l = get_param_from_dataset_name("loss_rate", dataset)
  t = get_param_from_dataset_name("transfer_rate", dataset)
  output = "../BenoitDatasets/families"
  jprime.generate_jprime(species_internal, families, sites, model, bl_factor, d, l, t, output, seed) 

def run_reference_methods(dataset):
  print("*************************************")
  print("Run reference methods for " + dataset)
  print("*************************************")
  dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
  is_dna = 1
  starting_trees = 10
  bs_trees = 100
  cores = 40
  save_sdtout = sys.stdout
  redirected_file = os.path.join(dataset_dir, "run_all.txt")
  print("Redirected logs to " + redirected_file)
  sys.stdout = open(redirected_file, "w")
  run_all.run_reference_methods(dataset_dir, is_dna, starting_trees, bs_trees, cores)
  sys.stdout = save_sdtout
  print("End of run_all")

def run_all_reference_methods(datasets):
  for dataset in datasets:
    run_reference_methods(dataset)

def run_all_generax(datasets, raxml = True, random = True):
  for dataset in datasets:
    print("*************************************")
    print("Run generate dataset for " + dataset)
    print("*************************************")
    is_dna = (not dataset in protein_datasets)
    if (raxml):
      run_generax(dataset, "raxml", False, "GeneRax-DL-Raxml", is_dna)
      run_generax(dataset, "raxml", True, "GeneRax-DTL-Raxml", is_dna)
    if (random):
      run_generax(dataset, "random", False, "GeneRax-DL-Random", is_dna)
      run_generax(dataset, "random", True, "GeneRax-DTL-Random", is_dna)

def generate_all_datasets(datasets):
  for dataset in datasets:
    generate_dataset(dataset)

def get_available_datasets(prefix):
  res = []
  for dataset in os.listdir(exp.families_datasets_root):
    try:
      if (not dataset.startswith(prefix)):
        continue
      res.append(dataset)
    except:
      continue
  return res

def get_rf_from_logs(logs):
  lines = logs.split("\n")
  dico = {}
  for line in lines:
    if (line.startswith("- True - ")):
      line = line.replace("\t", " ")
      split = line.split(" ")
      method = split[3][:-1]
      rf = split[-1]
      dico[method] = rf
  return dico


def get_timings(dataset):
  dico = {}
  lines = open(os.path.join(dataset, "runtimes.txt")).readlines()
  for line in lines:
    split = line.replace("\n", "").split(" ")
    dico[split[0]] = split[1]
  return dico

def get_results(dataset):
  try:
    analyse_script = os.path.join(exp.tools_root, "families", "analyze_dataset.py")
    families_path = os.path.join(exp.benoit_datasets_root, "families", dataset, "families")
    cmd = []
    cmd.append("python")
    cmd.append(analyse_script)
    cmd.append(families_path)
    logs = subprocess.check_output(cmd)
    return get_rf_from_logs(logs) 
  except:
    print("Failed to get RF distances from dataset " + dataset)
    return None



def get_datasets_to_plot(datasets_rf_dico, fixed_params_dico, x_param):

  datasets_to_plot = []
  for dataset in datasets_rf_dico:
    ok = True
    for fixed_param in fixed_params_dico:
      if (fixed_param == x_param):
        continue
      dataset_param_value = get_param_from_dataset_name(fixed_param, dataset)
      if (dataset_param_value != fixed_params_dico[fixed_param]):
        print("No: " + fixed_param + " " + dataset_param_value)
        ok = False
        continue
    if (ok):
      datasets_to_plot.append(dataset)
  return datasets_to_plot



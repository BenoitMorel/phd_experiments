import subprocess
import sys
import os 

sys.path.insert(0, os.path.join("tools", "families"))
import generate_families_with_jprime as jprime
import run_all

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
    try:
      return dataset.split("_")[8][1:]
    except:
      return "0.0"
  elif (parameter == "dl_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    res =  str(float(d)/float(l))
    print(dataset + " " + res)
    return res
  else:
    return "invalid"

# couples of (species interval, seed) to get a given number of species node
jsim_species_to_params = {}
jsim_species_to_params[25] = (8, 42)

def run_generax(dataset, starting_tree, with_transfers, run_name):
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
  print("Running " + " ".join(command))
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
  dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
  is_dna = 1
  starting_trees = 10
  bs_trees = 100
  cores = 40
  run_all.run_reference_methods(dataset_dir, is_dna, starting_trees, bs_trees, cores)

def run_all_reference_methods(datasets):
  for dataset in datasets:
    run_reference_methods(dataset)

def run_all_generax(datasets):
  for dataset in datasets:
    run_generax(dataset, "raxml", False, "GeneRax-DL-Raxml")
    run_generax(dataset, "random", False, "GeneRax-DL-Random")
    run_generax(dataset, "raxml", True, "GeneRax-DTL-Raxml")
    run_generax(dataset, "random", True, "GeneRax-DTL-Random")

def generate_all_datasets(datasets):
  for dataset in datasets:
    generate_dataset(dataset)



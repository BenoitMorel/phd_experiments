import re
import fam
import sys
import os
import generate_families_with_jprime as jprime
import generate_families_with_zombi as zombi
import generate_families_with_simphy as simphy
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp

# couples of (species interval, seed) to get a given number of species node
jsim_species_to_params = {}
jsim_species_to_params[5] = (3, 2)
jsim_species_to_params[10] = (3, 13)
jsim_species_to_params[12] = (3, 9)
jsim_species_to_params[14] = (3, 17)
jsim_species_to_params[19] = (3, 11)
jsim_species_to_params[27] = (3, 12)
jsim_species_to_params[41] = (3, 4)
jsim_species_to_params[16] = (3, 6)


def get_param_from_dataset_name(parameter, dataset):
  if (parameter == "tag"):
    return dataset.split("_")[1]
  elif (parameter == "species"):
    return dataset.split("_")[2][1:]
  elif (parameter == "families"):
    return dataset.split("_")[3][1:]
  elif (parameter == "sites"):
    return dataset.split("_")[4][5:]
  elif (parameter == "model"):
    return dataset.split("_")[5]
  elif (parameter == "bl"):
    return dataset.split("_")[6][2:]
  elif (parameter == "dup_rate"):
    return dataset.split("_")[7][1:]
  elif (parameter == "loss_rate"):
    return dataset.split("_")[8][1:]
  elif (parameter == "transfer_rate"):
    return dataset.split("_")[9][1:]
  elif (parameter == "gene_conversion_rate"):
    return dataset.split("_")[10][2:]
  elif (parameter == "perturbation"):
    return dataset.split("_")[11][1:]
  elif (parameter == "population"):
    return dataset.split("_")[12][3:]
  elif (parameter == "ms"):
    return dataset.split("_")[13][2:]
  elif (parameter == "msmf"):
    return dataset.split("_")[13][2:]
  elif (parameter == "mf"):
    return dataset.split("_")[14][2:]
  elif (parameter == "av_miss"):
    ms = float(get_param_from_dataset_name("ms", dataset))
    mf = float(get_param_from_dataset_name("mf", dataset))
    return str(round(ms * mf, 2))
  elif (parameter == "seed"):
    return dataset.split("_")[15][4:]
  elif (parameter == "tl_ratio"):
    t = get_param_from_dataset_name("transfer_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    if (float(l) == 0.0):
      return "-1.0"
    return  str(float(t)/float(l))
  elif (parameter == "dl_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    if (float(l) == 0.0):
      return "-1.0"
    return str(float(d)/float(l))
  elif (parameter == "dt_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    t = get_param_from_dataset_name("transfer_rate", dataset)
    if (float(t) == 0.0):
      return "-1.0"
    return str(float(d)/float(t))
  elif (parameter == "av_rate"):
    d = float(get_param_from_dataset_name("dup_rate", dataset))
    l = float(get_param_from_dataset_name("loss_rate", dataset))
    t = float(get_param_from_dataset_name("transfer_rate", dataset))
    if (float(t) == 0.0):
      return "-1.0"
    return str((d + t + l) / 2.0)
  elif (parameter == "discordance"):
    return float(fam.get_discordance_rate(fam.get_datadir(dataset)))
  else:
    return "invalid"

def generate_dataset(dataset):
  tag = get_param_from_dataset_name("tag", dataset)
  species = get_param_from_dataset_name("species", dataset)
  families = int(get_param_from_dataset_name("families", dataset))
  sites = get_param_from_dataset_name("sites", dataset)
  model = get_param_from_dataset_name("model", dataset)
  bl_factor = float(get_param_from_dataset_name("bl", dataset))
  d = get_param_from_dataset_name("dup_rate", dataset)
  l = get_param_from_dataset_name("loss_rate", dataset)
  t = get_param_from_dataset_name("transfer_rate", dataset)
  gc = get_param_from_dataset_name("gene_conversion_rate", dataset)
  p = get_param_from_dataset_name("perturbation", dataset)
   
  output = exp.families_datasets_root
  if (dataset.startswith("jsim")):
    species_internal, seed = jsim_species_to_params[int(species)]
    jprime.generate_jprime(tag, species_internal, families, sites, model, bl_factor, d, l, t, p, output, seed) 
  elif (dataset.startswith("zsim")):
    zombi.generate_zombi(tag, species, families, sites, model, bl_factor, d, l, t, output) 
  elif (dataset.startswith("ssim")):
    seed = get_param_from_dataset_name("seed", dataset)
    population = get_param_from_dataset_name("population", dataset)
    miss_species = get_param_from_dataset_name("ms", dataset)
    miss_fam = get_param_from_dataset_name("mf", dataset)
    model = "GTR"
    simphy.generate_simphy(tag, species, families, sites, model, bl_factor, d, l, t, gc, p, population, miss_species, miss_fam, output,  seed) 
  else:
    print("Unknown simulator for dataset " + dataset)
    sys.exit(1)

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

def get_param_position(fixed_point, param_name):
  split = fixed_point.split("_")
  if (param_name == "tag"):
    return 1
  for i in range(0, len(split)):
    if (param_name ==  re.sub("[0-9]*[\.]*[0-9]*", "", split[i])):
      return i
  print("ERROR: unknown parameter " + param_name)
  raise

def change_param_in_dataset_name(dataset, param_name, new_value):
  split = dataset.split("_")
  split[get_param_position(dataset, param_name)] = param_name + str(new_value)

def is_float(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

def get_dataset_variations(datasets, fixed_point, strings_to_replace):
  for string_to_replace in strings_to_replace:
    if (string_to_replace == "none"):
      datasets.append(fixed_point)
      continue
    elems_to_replace = string_to_replace.split("_")
    split = fixed_point.split("_")
    for elem in elems_to_replace:
      param_name =  re.sub("[0-9]*[\.]*[0-9]*", "", elem)
      param_name = re.sub("e-", "", param_name)
      param_value =  elem[len(param_name):]
      param_position = get_param_position(fixed_point, param_name)
      split[param_position] = param_name + str(param_value)
    dataset = "_".join(split)
    if (dataset in datasets):
      exit(1)
    datasets.append(dataset)

def duplicate_families(input_datadir, output_datadir, family):
  src = fam.get_family_path(input_datadir, family)
  dest = fam.get_family_path(output_datadir, family)
  shutil.copytree(src, dest)
  
def duplicate_families_symlink(input_datadir, output_datadir, family):
  fam.init_family_directories(output_datadir, family)
  exp.relative_symlink(fam.get_alignment(input_datadir, family), fam.get_alignment(output_datadir, family))
  exp.relative_symlink(fam.get_mappings(input_datadir, family), fam.get_mappings(output_datadir, family))
  input_gene_trees_dir = fam.get_gene_tree_dir(input_datadir, family)
  output_gene_trees_dir = fam.get_gene_tree_dir(output_datadir, family)
  for f in os.listdir(input_gene_trees_dir):
    if (f.startswith("true") or f.startswith("raxml") or "mrbayes" in f):
      exp.relative_symlink(os.path.join(input_gene_trees_dir, f), os.path.join(output_gene_trees_dir, f))
  



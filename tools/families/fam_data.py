import fam
import generate_families_with_jprime as jprime
import generate_families_with_zombi as zombi
import generate_families_with_simphy as simphy

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

#protein_datasets = ["swiss", "cyano_simulated", "cyano_empirical", "sub_t0.05_s0.5_cyano_empirical", "sub_t0.01_s0.2_ensembl_8880_15"]

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
  elif (parameter == "perturbation"):
    return dataset.split("_")[10][1:]
  elif (parameter == "population"):
    return dataset.split("_")[11][3:]
  elif (parameter == "sample_mu"):
    return dataset.split("_")[12][2:]
  elif (parameter == "sample_theta"):
    return dataset.split("_")[13][5:]
  elif (parameter == "seed"):
    return dataset.split("_")[14][4:]
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
    if (float(l) == 0.0):
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
    return float(fam.get_discordance_rate(get_datadir(dataset)))
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
  p = get_param_from_dataset_name("perturbation", dataset)
   
  output = "../BenoitDatasets/families"
  if (dataset.startswith("jsim")):
    species_internal, seed = jsim_species_to_params[int(species)]
    jprime.generate_jprime(species_internal, families, sites, model, bl_factor, d, l, t, p, output, seed) 
  elif (dataset.startswith("zsim")):
    zombi.generate_zombi(species, families, sites, model, bl_factor, d, l, t, output) 
  elif (dataset.startswith("ssim")):
    seed = get_param_from_dataset_name("seed", dataset)
    population = get_param_from_dataset_name("population", dataset)
    mu = get_param_from_dataset_name("sample_mu", dataset)
    theta = get_param_from_dataset_name("sample_theta", dataset)
    model = "GTR"
    simphy.generate_simphy(tag, species, families, sites, model, bl_factor, d, l, t, p, population, mu, theta, output,  seed) 
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
  exit(1)


import subprocess
import sys
import os 
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, 'scripts')
import experiments as exp
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import generate_families_with_jprime as jprime
import run_raxml_supportvalues as raxml
import run_ALE
import run_all
import run_phyldog_light
import saved_metrics
import eval_generax_likelihood

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



protein_datasets = ["swiss", "cyano_simulated", "sub_t0.05_s0.5_cyano_empirical", "sub_t0.01_s0.2_ensembl_8880_15"]


def run_generax(dataset, starting_tree, with_transfers, run_name, is_dna, cores = 40):
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
  command.append("--run")
  command.append(run_name)
  if (not is_dna):
    command.append("--protein")
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)


def generate_dataset(dataset):
  species = fam.get_param_from_dataset_name("species", dataset)
  species_internal, seed = jsim_species_to_params[int(species)]
  families = int(fam.get_param_from_dataset_name("families", dataset))
  sites = fam.get_param_from_dataset_name("sites", dataset)
  model = fam.get_param_from_dataset_name("model", dataset)
  bl_factor = float(fam.get_param_from_dataset_name("bl", dataset))
  d = fam.get_param_from_dataset_name("dup_rate", dataset)
  l = fam.get_param_from_dataset_name("loss_rate", dataset)
  t = fam.get_param_from_dataset_name("transfer_rate", dataset)
  output = "../BenoitDatasets/families"
  jprime.generate_jprime(species_internal, families, sites, model, bl_factor, d, l, t, output, seed) 

def run_reference_methods(dataset, cores = 40):
  print("*************************************")
  print("Run reference methods for " + dataset)
  print("*************************************")
  dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
  is_dna = (not dataset in protein_datasets)
  starting_trees = 10
  bs_trees = 100
  save_sdtout = sys.stdout
  redirected_file = os.path.join(dataset_dir, "logs_run_all.txt")
  print("Redirected logs to " + redirected_file)
  sys.stdout.flush()
  sys.stdout = open(redirected_file, "w")
  run_all.run_reference_methods(dataset_dir, is_dna, starting_trees, bs_trees, cores)
  sys.stdout = save_sdtout
  print("End of run_all")
  sys.stdout.flush()

def run_all_phyldog(datasets, is_dna, cores):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    is_dna = (not dataset in protein_datasets)
    run_phyldog_light.run_phyldog_on_families(dataset_dir, is_dna, cores)

def run_all_raxml_light(datasets, cores = 40):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    is_dna = (not dataset in protein_datasets)
    starting_trees = 1
    bs_trees = 0
    raxml.run_pargenes_and_extract_trees(dataset_dir, is_dna, starting_trees, bs_trees, cores, "RAxML-light", False)

def run_all_ALE(datasets, is_dna, cores = 40):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_ALE.run_exabayes_and_ALE(dataset_dir, is_dna, cores)

def run_all_reference_methods(datasets, cores = 40):
  for dataset in datasets:
    run_reference_methods(dataset, cores)

def run_all_generax(datasets, raxml = True, random = True, DL = True, DTL = True, cores = 40):
  for dataset in datasets:
    print("*************************************")
    print("Run generate dataset for " + dataset)
    print("*************************************")
    is_dna = (not dataset in protein_datasets)
    if (raxml):
      if (DL):
        run_generax(dataset, "RAxML-NG", False, "GeneRax-DL-Raxml", is_dna, cores)
      if (DTL):
        run_generax(dataset, "RAxML-NG", True, "GeneRax-DTL-Raxml", is_dna, cores)
    if (random):
      if (DL):
        run_generax(dataset, "Random", False, "GeneRax-DL-Random", is_dna, cores)
      if (DTL):
        run_generax(dataset, "Random", True, "GeneRax-DTL-Random", is_dna, cores)

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


def get_runtimes(dataset):
  dico = {}
  try:
    lines = open(os.path.join(exp.benoit_datasets_root, "families", dataset, "runtimes.txt")).readlines()
    for line in lines:
      split = line.replace("\n", "").split(" ")
      dico[split[0]] = split[1]
  except:
    return None
  return dico

def get_results(dataset):
  #try:
    analyse_script = os.path.join(exp.tools_root, "families", "analyze_dataset.py")
    dataset_path = os.path.join(exp.benoit_datasets_root, "families", dataset)
    cmd = []
    cmd.append("python")
    cmd.append(analyse_script)
    cmd.append(dataset_path)
    print(" ".join(cmd))
    logs = subprocess.check_output(cmd).decode("utf-8")
    return get_rf_from_logs(logs) 
  #except:
  #  print("Failed to get RF distances from dataset " + dataset)
  #  return None

def run_all_analyzes(datasets):
  for dataset in datasets:
    print("analyze " + dataset)
    get_results(dataset)

def compute_likelihoods(datasets, cores = 40):
  for dataset in datasets:
    dataset_dir = os.path.join(exp.benoit_datasets_root, "families", dataset)
    is_protein = dataset in protein_datasets
    starting_trees = ["RAxML-NG", "Treerecs", "ALE-DL", "ALE-DTL", "GeneRax-DL-Random", "GeneRax-DTL-Random", "True", "Notung"]
    #starting_trees = ["GeneRax-DL-Raxml", "GeneRax-DTL-Raxml"]
    for tree in starting_trees:
      eval_generax_likelihood.eval_and_save_likelihood(dataset_dir, tree, 0, is_protein, cores)  
      eval_generax_likelihood.eval_and_save_likelihood(dataset_dir, tree, 1, is_protein, cores)  

def get_datasets_to_plot(datasets_rf_dico, fixed_params_dico, x_param):

  datasets_to_plot = []
  for dataset in datasets_rf_dico:
    ok = True
    for fixed_param in fixed_params_dico:
      if (fixed_param == x_param):
        continue
      dataset_param_value = fam.get_param_from_dataset_name(fixed_param, dataset)
      if (dataset_param_value != fixed_params_dico[fixed_param]):
        ok = False
        continue
    if (ok):
      datasets_to_plot.append(dataset)
  return datasets_to_plot

def get_metrics_for_datasets(datasets_prefix, metric_name):
  datasets = get_available_datasets(datasets_prefix)
  datasets_rf_dico = {}
  datasets_runtimes_dico = {}
  total = len(datasets)
  for dataset in datasets:
    #if (metric_name == "average_rrf"):
    #  get_results(dataset)
    dataset_dir = os.path.join(exp.families_datasets_root, dataset)
    res = saved_metrics.get_metrics(dataset_dir, metric_name)
    if (res != None):
      if (metric_name == "runtimes"):
        for key in res: 
          if ("ALE" in key):
            res[key] = str(float(res[key]) + float(res["ExaBayes"]))
      datasets_rf_dico[dataset] = res
  return datasets_rf_dico


def get_param_position(fixed_point, param_name):
  split = fixed_point.split("_")
  for i in range(0, len(split)):
    if (param_name ==  re.sub("[0-9]*[\.]*[0-9]*", "", split[i])):
      return i
  print("ERROR: unknown parameter " + param_name)
  exit(1)

def add_dataset(datasets, fixed_point, strings_to_replace):
  for string_to_replace in strings_to_replace:
    elems_to_replace = string_to_replace.split("_")
    split = fixed_point.split("_")
    for elem in elems_to_replace:
      param_name =  re.sub("[0-9]*[\.]*[0-9]*", "", elem)
      param_value =  re.sub("[a-zA-Z]*", "", elem)
      param_position = get_param_position(fixed_point, param_name)
      split[param_position] = param_name + str(param_value)
    dataset = "_".join(split)
    print("Add " + dataset)
    if (dataset in datasets):
      print("duplicate: " + dataset)
      exit(1)
    datasets.append(dataset)


def plot(datasets_rf_dico, x_param, fixed_params_dico, methods, x_label, y_label, output):
  print(x_param)
  datasets_to_plot = get_datasets_to_plot(datasets_rf_dico, fixed_params_dico, x_param)
  #print(datasets_to_plot)
  datasets_to_plot.sort(key = lambda t: float(fam.get_param_from_dataset_name(x_param, t)))
  df = pd.DataFrame()
  f, ax = plt.subplots(1)
  fake_df = {}
  fake_df[x_param] = []
  for method in methods:
    fake_df[method] = []
  for dataset in datasets_to_plot:
    #print(float(fam.get_param_from_dataset_name(x_param,     dataset)))
    fake_df[x_param].append(float(fam.get_param_from_dataset_name(x_param, dataset)))
    rf_dico = datasets_rf_dico[dataset]
    for method in methods:
      if (not method in rf_dico):
        print("Warning: missing data for method " + method + " and dataset " + dataset)
    for method in rf_dico:
      if (not method in methods):
        #print("Unknown method " + method)
        continue
      fake_df[method].append(float(rf_dico[method]))
  for elem in fake_df:
    print("\t" + str(elem))
    print("\t" + str(fake_df[elem]))
    df[elem] = fake_df[elem]
  
  for method in methods:
    style = "solid"
    plt.plot(x_param, method, data=df, marker='x', linestyle = style, linewidth=2, label = method)
    #ax.set_ylim(bottom=0)
  plt.xlabel(x_label[x_param])
  plt.ylabel(y_label)
  plt.legend()
  plt.savefig(output)
  print("Saving result in " + output)
  plt.close()


class Plotter(object):
  def __init__(self, datasets_values_dico, fixed_parameters, methods, x_labels, y_label, prefix):
    self.datasets_values_dico = datasets_values_dico
    self.fixed_parameters = fixed_parameters
    self.methods = methods
    self.x_labels = x_labels
    self.prefix = prefix
    self.y_label = y_label
  
  def __call__(self, parameter):
    plot(self.datasets_values_dico, parameter, self.fixed_parameters, self.methods, self.x_labels,  self.y_label, self.prefix + "_" + parameter + ".png")

def submit_single_experiment_haswell(dataset, do_generate, cores):
  command = []
  command.append("python")
  command.append(os.path.join(exp.tools_root, "publications_generax", "single_experiment.py"))
  command.append(dataset)
  command.append(str(do_generate))
  command.append(str(cores))
  results_dir = os.path.join("single_experiments", dataset)
  results_dir = exp.create_result_dir(results_dir, [])
  result_msg = ""
  exp.write_results_info(results_dir, result_msg) 
  submit_path = os.path.join(results_dir, "sub_generax.sh")
  exp.submit(submit_path, " ".join(command), cores, "haswell") 

def submit_multiple_experiments_haswell(datasets, do_generate, cores):
  for dataset in datasets:
    submit_single_experiment_haswell(dataset, do_generate, cores)
  



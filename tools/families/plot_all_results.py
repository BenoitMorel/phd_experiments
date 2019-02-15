import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp

possible_species = ["8"]
possible_sites = ["250", "500", "1000", "1500"] 
possible_bl = ["0.5", "1.0", "2.0", "4.0"]
possible_dup_rates = ["1.0", "0.5", "0.1"]

possible_parameters = {}
possible_parameters["species"] = possible_species 
possible_parameters["sites"] = possible_sites 
possible_parameters["bl"] = possible_bl
possible_parameters["dup_rates"] = possible_dup_rates


def parameters_to_dataset(species, sites, bl, dup_rate):
  loss_rate = str(float(rates) / 2.0)
  return "jsim_s" + species + "_f50_sites" + sites + "_dna4_bl" + bl + "_d" + dup_rate + "_" + loss_rate + "_seed42"

def get_param_from_dataset_name(parameter, dataset):
  if (parameter == "species"):
    return dataset.split("_")[1][1:]
  elif (parameter == "sites"):
    return dataset.split("_")[3][5:]
  elif (parameter == "bl"):
    return dataset.split("_")[5][2:]
  elif (parameter == "dup_rates"):
    return dataset.split("_")[6][1:]
  else:
    return "invalid"

def get_available_datasets(prefix):
  res = []
  for dataset in os.listdir(exp.families_datasets_root):
    try:
      ok = True
      if (not dataset.startswith(prefix)):
        continue
      for parameter in possible_parameters:
        if (not get_param_from_dataset_name(parameter, dataset) in possible_parameters[parameter]):
          ok = False
          continue
      if (ok):
        res.append(dataset)
      pass
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

def get_results(dataset):
  analyse_script = os.path.join(exp.tools_root, "families", "analyse_dataset.py")
  families_path = os.path.join(exp.benoit_datasets_root, "families", dataset, "families")
  results_path = os.path.join(exp.results_root, "MultipleJointSearch", dataset, "SPR_2_start_raxml_split", "normald_40", "run_0", "scheduler_run")
  cmd = []
  cmd.append("python")
  cmd.append(analyse_script)
  cmd.append(families_path)
  cmd.append(results_path)
  logs = subprocess.check_output(cmd)
  return get_rf_from_logs(logs) 


def get_datasets_to_plot(datasets_rf_dico, fixed_params_dico):

  datasets_to_plot = []
  for dataset in datasets_rf_dico:
    ok = True
    for fixed_param in fixed_params_dico:
      dataset_param_value = get_param_from_dataset_name(fixed_param, dataset)
      if (dataset_param_value != fixed_params_dico[fixed_param]):
        ok = False
        continue
    if (ok):
      datasets_to_plot.append(dataset)
  return datasets_to_plot


def plot(datasets_rf_dico, x_param, fixed_params_dico, output):
  datasets_to_plot = get_datasets_to_plot(datasets_rf_dico, fixed_params_dico)
  datasets_to_plot.sort(key = lambda t: float(get_param_from_dataset_name(x_param, t)))
  df = pd.DataFrame()
  f, ax = plt.subplots(1)
  ax.set_ylim(ymin=0)
  
  methods = []
  fake_df = {}
  fake_df[x_param] = []
  for dataset in datasets_to_plot:
    fake_df[x_param].append(float(get_param_from_dataset_name(x_param, dataset)))
    rf_dico = datasets_rf_dico[dataset]
    for method in rf_dico:
      if (not method in methods):
        methods.append(method)
        fake_df[method] = []
      fake_df[method].append(float(rf_dico[method]))
  for elem in fake_df:
    df[elem] = fake_df[elem]
  
  for method in methods:
    plt.plot(x_param, method, data=df, marker='x', linewidth=2)

  plt.xlabel(x_param)
  plt.ylabel('RF distance')
  plt.legend()
  plt.savefig(output)
  plt.close()

datasets = get_available_datasets("jsim")
datasets_rf_dico = {}
for dataset in datasets:
  datasets_rf_dico[dataset] = get_results(dataset)



params_value_dico_sites = {}
params_value_dico_sites["species"] = "8"
params_value_dico_sites["dup_rates"] = "1.0"
params_value_dico_sites["bl"] = "1.0"
plot(datasets_rf_dico, "sites", params_value_dico_sites, "sites.png")

params_value_dico_sites = {}
params_value_dico_sites["species"] = "8"
params_value_dico_sites["sites"] = "500"
params_value_dico_sites["bl"] = "1.0"
plot(datasets_rf_dico, "dup_rates", params_value_dico_sites, "rates.png")

params_value_dico_sites = {}
params_value_dico_sites["species"] = "8"
params_value_dico_sites["sites"] = "500"
params_value_dico_sites["dup_rates"] = "0.5"
plot(datasets_rf_dico, "bl", params_value_dico_sites, "bl.png")






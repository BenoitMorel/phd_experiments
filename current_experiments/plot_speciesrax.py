import os
import sys
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import experiments as exp
import fam
import fam_data
import saved_metrics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style("darkgrid")




def get_df(df):
  res = pd.DataFrame()
  for elem in df:
    res[elem] = df[elem]
  return res

def get_default_value(metric_name):
  if ("time" in metric_name):
    return 0.0
  if ("rf" in metric_name):
    return 1.0
  return -1.0

def plot(grouped_datasets, x_param, methods, subst_model, metric_name, output):
  df = {}
  datasets_keys = []
  method_dict = {}
  for method in methods:
    method_dict[method] = method
  method_dict["speciesrax-dtl-raxml-HYBRID"] = "generax"
  for key in grouped_datasets:
    datasets_keys.append(key)
  print("Sorting by " + x_param)
  datasets_keys.sort(key = lambda t: float(fam_data.get_param_from_dataset_name(x_param, t)))
  print(datasets_keys)
  f, ax = plt.subplots(1)
  df[x_param] = []
  for method in methods:
    df[method] = []
  for dataset_key in datasets_keys:
    # value for the varying parameter
    df[x_param].append(fam_data.get_param_from_dataset_name(x_param, grouped_datasets[dataset_key][0]))
    for method in methods:
      key = (method + "." + subst_model).lower()
      average = 0.0
      for dataset in grouped_datasets[dataset_key]:
        dataset_dir = os.path.join(exp.families_datasets_root, dataset)
        metrics = saved_metrics.get_metrics(dataset_dir, metric_name)
        if (metrics != None and key in metrics):
          average += float(metrics[key])
        else:
          average += float(get_default_value(metric_name))
          print("WARNING: missing value for run " + key + " and dataset " + dataset)
      average /= float(len(grouped_datasets[dataset_key]))
      df[method].append(average)
  df = get_df(df)
  for method in methods:
    print(x_param)
    print(method)
    print(df)
    #plt.plot(x_param, method, data=df, marker='.', linestyle = "solid", linewidth=2, label = method, markersize=12)
    plt.plot(x_param, method, data=df, marker='.', linewidth=2, label = method_dict[method], markersize=12)
  plt.xlabel(x_param)
  plt.ylabel(metric_name)
  plt.legend()
  plt.savefig(output)
  print("Saving result in " + output)
  plt.close()



def get_relevant_datasets(datasets, param, fixed_params_values):
  res = []
  for dataset in datasets:
    ok = True
    for fixed_param in fixed_params_values:
      if (fixed_param == param):
        continue
      value = fam_data.get_param_from_dataset_name(fixed_param, dataset)
      if (value != fixed_params_values[fixed_param]):
        print (str(fixed_param) + " " + str(value) + " " + str(fixed_params_values[fixed_param]))
        ok = False
        break
    if (ok):
      res.append(dataset)
  return res

def merge_datasets_per_seed(datasets):
  grouped_datasets = {}
  print("yop + " + str(datasets))
  for dataset in datasets:
    key = fam.get_first_dataset_starting_with(dataset.split("seed")[0])
    #key = fam.get_datadir(dataset)
    print(key)
    if (not key in grouped_datasets):
      grouped_datasets[key] = []
    grouped_datasets[key].append(dataset)
  return grouped_datasets

def get_plot_name(simulation_name, param, subst_model, metric_name):
  return "plot_" + simulation_name + "_" + metric_name.replace("_", "-") + "_" + param + "_" + subst_model

def plot_metric(param, fixed_params_values, methods, subst_model, metric_name, datasets, plot_name):
  relevant_datasets = get_relevant_datasets(datasets, param, fixed_params_values)
  grouped_datasets = merge_datasets_per_seed(relevant_datasets)
  plot(grouped_datasets, param, methods, subst_model, metric_name, plot_name + ".svg")
  


def get_datasets(prefix):
  all_datasets = os.listdir(os.path.join(exp.benoit_datasets_root, "families"))
  res = []
  for dataset in all_datasets:
    if (dataset.startswith(prefix)):
      res.append(dataset)
  return res
  
def generate_plot(datasets, params_to_plot, metric_names, methods, simulation_name, fixed_params_values, subst_model):
  for param in params_to_plot:
    for metric_name in metric_names:
      plot_name = get_plot_name(simulation_name, param, subst_model, metric_name)
      plot_metric(param, fixed_params_values, methods, subst_model, metric_name, datasets, plot_name)

def plot_parameters(methods, metric_names, param):
  simulation_name = "params"
  subst_model = "GTR"
  datasets = get_datasets("ssim")
  fixed_params_values = {}
  fixed_params_values["species"] = "40"
  fixed_params_values["bl"] = "1.0"
  fixed_params_values["families"] = "100"
  fixed_params_values["sites"] = "100"
  if (param == "transfer_rate"):
    fixed_params_values["tl_ratio"] = "1.0"
  if (param == "dup_rate" or param == "discordance"):
    fixed_params_values["transfer_rate"] = "0.0"
  else:
    fixed_params_values["transfer_rate"] = "0.1"
  if (param == "discordance"):
    fixed_params_values["dup_rate"] = "0.0"

  if (param == "transfer_rate" or param == "dup_rate"):
    fixed_params_values["population"] = "10"
  elif (param != "discordance"):
    fixed_params_values["population"] = "10000"
  generate_plot(datasets, [param], metric_names, methods, simulation_name, fixed_params_values, subst_model)
  
def plot_missing(methods, metric_names):
  simulation_name = "missing"
  subst_model = "GTR"
  datasets = get_datasets("ssim")
  fixed_params_values = {}
  fixed_params_values["species"] = "20"
  fixed_params_values["bl"] = "1.0"
  fixed_params_values["families"] = "100"
  fixed_params_values["sites"] = "100"
  fixed_params_values["transfer_rate"] = "0.1"
  fixed_params_values["dup_rate"] = "0.1"
  fixed_params_values["population"] = "10000"
  fixed_params_values["sample_mu"] = "0.7"
  fixed_params_values["sample_theta"] = "0.0"
  generate_plot(datasets, ["sample_mu"], metric_names, methods, simulation_name, fixed_params_values, subst_model)
  generate_plot(datasets, ["sample_theta"], metric_names, methods, simulation_name, fixed_params_values, subst_model)
  fixed_params_values["sample_mu"] = "0.35"
  simulation_name = "missingFamMu"
  generate_plot(datasets, ["families"], metric_names, methods, simulation_name, fixed_params_values, subst_model)
  simulation_name = "missingFamTheta"
  fixed_params_values["sample_mu"] = "0.75"
  fixed_params_values["sample_theta"] = "0.15"
  generate_plot(datasets, ["families"], metric_names, methods, simulation_name, fixed_params_values, subst_model)
  

def plot_species(methods, metric_names):
 
  #only one can be True:
  subst_model = "GTR"
  dtl_mode = False
  dl_mode = False
  idtl_mode = False
  i_mode = False
  
  dtl_mode = True
  simulation_name = ""
  if (dtl_mode):
    simulation_name += "dtl"
  elif (dl_mode):
    simulation_name += "dl"
  elif (idtl_mode):
    simulation_name += "idtl"
  elif (i_mode):
    simulation_name += "ils"
  datasets = get_datasets("ssim") 
  #params_to_plot = ["dup_rate"]
  params_to_plot = ["species"]
  fixed_params_values_dtl = {}
  fixed_params_values_dtl["species"] = "40"
  fixed_params_values_dtl["bl"] = "1.0"
  fixed_params_values_dtl["families"] = "100"
  fixed_params_values_dtl["sites"] = "100"
  fixed_params_values_dtl["tag"] = simulation_name
  if (dtl_mode or dl_mode):
    fixed_params_values_dtl["dup_rate"] = "0.2"
  elif (idtl_mode):
    fixed_params_values_dtl["dup_rate"] = "0.1"
  else:
    fixed_params_values_dtl["dup_rate"] = "0.0"
  if (dtl_mode):
    fixed_params_values_dtl["transfer_rate"] = "0.2"
  elif (idtl_mode):
    fixed_params_values_dtl["transfer_rate"] = "0.1"
  else:
    fixed_params_values_dtl["transfer_rate"] = "0.0"
  if (i_mode):
    fixed_params_values_dtl["population"] = "100000"
  elif (idtl_mode):
    fixed_params_values_dtl["population"] = "20000"
  else:
    fixed_params_values_dtl["population"] = "10"
  generate_plot(datasets, params_to_plot, metric_names, methods, simulation_name, fixed_params_values_dtl, subst_model)

def main_plot_metrics():
  #methods = ["duptree", "astral", "astralpro", "njrax-NJst", "speciesrax-dtl-raxml-HYBRID", "concatenation-naive"]
  methods = ["astral", "astralpro", "speciesrax-dtl-raxml-HYBRID"]
  #methods = ["duptree", "njst-original", "njst-reweighted", "astralpro", "njrax-NJst", "speciesrax-dtl-raxml-HYBRID", "speciesrax-dtl-raxml-SPR", "speciesrax-prune-SPR", "speciesrax-per-family-SPR"]
  metric_names = ["species_unrooted_rf",  "runtimes"]
  plot_species(methods, metric_names)
  #plot_parameters(methods, metric_names, "dup_rate")
  #plot_parameters(methods, metric_names, "transfer_rate")
  #plot_parameters(methods, metric_names, "bl")
  #plot_parameters(methods, metric_names, "families")
  #plot_parameters(methods, metric_names, "discordance")
  #metric_names = ["species_unrooted_rf"]
  #plot_missing(methods, metric_names)

if (__name__ == "__main__"):
  main_plot_metrics()


import os
import sys
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/plotters')
import experiments as exp
import fam
import fam_data
import saved_metrics

import plot_histogram
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

"""
  grouped_datasets: dict, with key being a dataset name, and the value a list of this dataset with different seeds
  methods: list of methods to plot
  methods_dict: dictionary that mapped a method to the name that should be displaed
  subst_model:
  metric_name: metric to plot (for instance runtimes)
  output: output file
"""
def plot(grouped_datasets, x_param, methods, methods_dict, subst_model, metric_name, output, get_param_fct = fam_data.get_param_from_dataset_name, title = None):
  df = {}
  datasets_keys = []
  method_dict = {}
  for method in methods:
    if (not method in methods_dict):
      methods_dict[method] = method
  for key in grouped_datasets:
    datasets_keys.append(key)
  print("Sorting by " + x_param)
  print(datasets_keys)
  datasets_keys.sort(key = lambda t: float(get_param_fct(x_param, t)))
  print(datasets_keys)
  f, ax = plt.subplots(1)
  df[x_param] = []
  for method in methods:
    df[method] = []
  for dataset_key in datasets_keys:
    # value for the varying parameter
    df[x_param].append(get_param_fct(x_param, grouped_datasets[dataset_key][0]))
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
    #print(x_param)
    #print(method)
    #print(df)
    method_alias = methods_dict[method][0]
    linestyle = None
    color = None
    if (len(methods_dict[method]) > 1):
      linestyle = methods_dict[method][1]
    if (len(methods_dict[method]) > 2):
      color = methods_dict[method][2]
    
    method_marker = "."
    markersize = 12
    if (linestyle != None and "dash" in linestyle):
      method_marker = "x"
      markersize = 8
    #plt.plot(x_param, method, data=df, marker='.', linestyle = "solid", linewidth=2, label = method, markersize=12)
    plt.plot(x_param, method, data=df, marker=method_marker, linewidth=2, label = method_alias, markersize=markersize, linestyle = linestyle, color = color)
  plt.xlabel(x_param)
  plt.ylabel(metric_name)
  if (title != None):
    plt.title(title)
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
      if (fixed_param == "population" and param == "discordance"):
        continue
      value = fam_data.get_param_from_dataset_name(fixed_param, dataset)
      if (value != fixed_params_values[fixed_param]):
        #print("get param from " + fixed_param + " " + dataset + ": " + value)
        #print (str(fixed_param) + " " + str(value) + " " + str(fixed_params_values[fixed_param]))
        ok = False
        break
    if (ok):
      res.append(dataset)
  return res

def merge_datasets_per_seed(datasets):
  grouped_datasets = {}
  for dataset in datasets:
    key = fam.get_first_dataset_starting_with(dataset.split("seed")[0])
    #key = fam.get_datadir(dataset)
    #print(key)
    if (not key in grouped_datasets):
      grouped_datasets[key] = []
    grouped_datasets[key].append(dataset)
  return grouped_datasets

def get_plot_name(simulation_name, param, subst_model, metric_name):
  return "plot_" + simulation_name + "_" + metric_name.replace("_", "-") + "_" + param + "_" + subst_model

def plot_metric(param, fixed_params_values, methods, methods_dict, subst_model, metric_name, datasets, plot_name):
  relevant_datasets = get_relevant_datasets(datasets, param, fixed_params_values)
  grouped_datasets = merge_datasets_per_seed(relevant_datasets)
  plot(grouped_datasets, param, methods, methods_dict, subst_model, metric_name, plot_name + ".svg")
  
def plot_runtimes(dataset, subst_model, methods, methods_dict):
  dataset_dir = fam.get_datadir(dataset)
  metrics = saved_metrics.get_metrics(dataset_dir, "runtimes")
  xlabels = []
  yvalues = []
  for method in methods:
    metric_key = (method + "." + subst_model).lower()
    if (metrics != None and metric_key in metrics):
      xlabels.append(methods_dict[method])
      yvalues.append(float(metrics[metric_key]))
    else:
      print("Error: cannot get runtime from " + dataset + " and " + metric_key)
  output = "plot_runtime_" + dataset + "_" + subst_model + ".svg"
  ycaption = "Runtimes (s)" 
  plot_histogram.plot_histogram(xlabels, yvalues, ycaption = ycaption, log_scale = True, output = output)



def get_datasets(prefix):
  all_datasets = os.listdir(os.path.join(exp.benoit_datasets_root, "families"))
  res = []
  for dataset in all_datasets:
    if (dataset.startswith(prefix)):
      res.append(dataset)
  return res
  
def generate_plot(datasets, params_to_plot, metric_names, methods, methods_dict, simulation_name, fixed_params_values, subst_model):
  for param in params_to_plot:
    for metric_name in metric_names:
      plot_name = get_plot_name(simulation_name, param, subst_model, metric_name)
      plot_metric(param, fixed_params_values, methods, methods_dict, subst_model, metric_name, datasets, plot_name)

def plot_params(methods, metric_names):
  simulation_name = "params"
  subst_model = "GTR+G"
  datasets = get_datasets("ssim")
  params_to_plot = ["species", "sites", "dup_rate", "families", "tl_ratio", "discordance"]
  fixed_params_values = {}
  fixed_params_values["species"] = "25"
  fixed_params_values["bl"] = "1.0"
  fixed_params_values["families"] = "100"
  fixed_params_values["sites"] = "100"
  fixed_params_values["tag"] = "var"
  fixed_params_values["dup_rate"] = "1.0"
  fixed_params_values["tl_ratio"] = "1.0"
  fixed_params_values["population"] = "470000000"
  generate_plot(datasets, params_to_plot, metric_names, methods, simulation_name, fixed_params_values, subst_model)

def main_plot_metrics():
  #methods = ["astralpro", "njrax-NJst", "speciesrax-dtl-raxml-HYBRID", "generax-select", "speciesrax-dtl-raxml-perfam-HYBRID"]
  methods = ["astralpro", "njrax-NJst", "speciesrax-dtl-raxml-perfam-HYBRID", "orthogenerax-all-njrax-njst", "generax-select", "generax-select-true"]
  metric_names = ["species_unrooted_rf", "species_rooted_rf"]
  plot_params(methods, metric_names)
  metric_names = ["species_rooted_rf"]#,  "runtimes"]
  plot_params(methods, metric_names)

if (__name__ == "__main__"):
  main_plot_metrics()


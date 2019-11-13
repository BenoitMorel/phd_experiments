import os
import sys
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import experiments as exp
import fam
import saved_metrics
import common

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
  for key in grouped_datasets:
    datasets_keys.append(key)
  datasets_keys.sort(key = lambda t: float(fam.get_param_from_dataset_name(x_param, t)))
  f, ax = plt.subplots(1)
  df[x_param] = []
  for method in methods:
    df[method] = []
  for dataset_key in datasets_keys:
    # value for the varying parameter
    df[x_param].append(fam.get_param_from_dataset_name(x_param, grouped_datasets[dataset_key][0]))
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
    plt.plot(x_param, method, data=df, marker='.', linewidth=2, label = method, markersize=12)
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
      value = fam.get_param_from_dataset_name(fixed_param, dataset)
      if (value != fixed_params_values[fixed_param]):
        ok = False
        break
    if (ok):
      res.append(dataset)
  return res

def merge_datasets_per_seed(datasets):
  grouped_datasets = {}
  for dataset in datasets:
    key = dataset.split("seed")[0]
    if (not key in grouped_datasets):
      grouped_datasets[key] = []
    grouped_datasets[key].append(dataset)
  return grouped_datasets

def get_plot_name(param, subst_model, metric_name):
  return "plot_" + metric_name.replace("_", "-") + "_" + param + "_" + subst_model

def plot_metric(param, fixed_params_values, methods, subst_model, metric_name, datasets):
  relevant_datasets = get_relevant_datasets(datasets, param, fixed_params_values)
  grouped_datasets = merge_datasets_per_seed(relevant_datasets)
  plot_name = get_plot_name(param, subst_model, metric_name)
  print(plot_name)
  plot(grouped_datasets, param, methods, subst_model, metric_name, plot_name + ".png")
  


def get_datasets(prefix):
  all_datasets = os.listdir(os.path.join(exp.benoit_datasets_root, "families"))
  res = []
  for dataset in all_datasets:
    if (dataset.startswith(prefix)):
      res.append(dataset)
  return res

def main_plot_metrics():
  datasets = get_datasets("ssim") 
  params_to_plot = ["species"]
  fixed_params_values_dtl = {}
  fixed_params_values_dtl["species"] = "20"
  fixed_params_values_dtl["dup_rate"] = "0.2"
  fixed_params_values_dtl["loss_rate"] = "0.2"
  fixed_params_values_dtl["transfer_rate"] = "0.2"
  fixed_params_values_dtl["bl"] = "1.0"
  fixed_params_values_dtl["families"] = "100"
  fixed_params_values_dtl["sites"] = "100"
  fixed_params_values_dtl["population"] = "10"
  
  methods = ["astral", "duptree", "stag", "speciesrax-dtl-raxml-NJ"]
  #methods = ["phyldogspecies", "duptree", "stag", "astral", "speciesrax-dtl-raxml-NJ", "speciesrax-dtl-raxml-slow-NJ"]
  #methods = ["speciesrax-dtl-raxml-slow", "speciesrax-dtl-raxml-slow-NJ"]
  subst_model = "GTR"
  metric_names = ["species_unrooted_rf", "runtimes"]

  for param in params_to_plot:
    for metric_name in metric_names:
      plot_metric(param, fixed_params_values_dtl, methods, subst_model, metric_name, datasets)


if (__name__ == "__main__"):
  main_plot_metrics()


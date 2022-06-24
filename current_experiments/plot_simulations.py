import os
import sys
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import fam_data
import fam
import experiments as exp
import saved_metrics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style("darkgrid")

def _get_param_to_datasets(datasets, param_name):
  param_to_datasets_dict = {}
  for dataset in datasets:
    p = str(fam_data.get_param_from_dataset_name(param_name, dataset))
    if (not p in param_to_datasets_dict):
      param_to_datasets_dict[p] = []
    param_to_datasets_dict[p].append(dataset)
  param_to_datasets = []
  for key in param_to_datasets_dict:
    pair = (key, param_to_datasets_dict[key])
    print(key + ": " + str(len(param_to_datasets_dict[key])) + " datasets")
    param_to_datasets.append(pair) 
  return param_to_datasets

def _get_average_methods_values(param_name, datasets, metric, methods, subst_model):
  values = {}
  for method in methods:
    values[method] = 0.0
    key = (method + "." + subst_model).lower()
    average = 0.0
    for dataset in datasets:
      dataset_dir = os.path.join(exp.families_datasets_root, dataset)
      metrics = saved_metrics.get_metrics(dataset_dir, metric)
      if (metrics != None and key in metrics):
        average += float(metrics[key])
      else:
        print("WARNING: missing value for run " + key + " and dataset " + dataset)
        average += float(get_default_value(metric_name))
    average /= float(len(datasets))
    values[method] = average
  return values

# values[param][method] = value
def _get_df(values, methods, param_name):
  temp = {}
  temp[param_name] = []
  for method in methods:
    temp[method] = []
  keys = list(values)
  keys.sort(key = lambda t: float(t))
  for param in keys:
    temp[param_name].append(param)
    for method in values[param]:
      temp[method].append(float(values[param][method]))
  res = pd.DataFrame()
  for elem in temp:
    res[elem] = temp[elem]
  return res
  
  
def _get_methods(methods_tuples):
  res = []
  for t in methods_tuples:
    res.append(t[0])
  return res

"""
  Generate a line plot for a varying simulation parameter
  - datasets: list of all datasets to include
  - param_name: the parameter to vary
  - metric: the metric to plot
  - method_tuples: list of tuple (method_key, name_to_display)
  - output: plot filename
"""
def plot_varying_params(datasets, param_name, metric, method_tuples, subst_model, output, xlabel = None, ylabel = None):
  param_to_datasets = _get_param_to_datasets(datasets, param_name)
  methods = _get_methods(method_tuples)
  values = {}  
  for xvalue,xdatasets in param_to_datasets:
    values[xvalue] = _get_average_methods_values(param_name, xdatasets, metric, methods, subst_model)
  df = _get_df(values, methods, param_name)
  for method_tuple in method_tuples:
    method = method_tuple[0]
    method_alias = method_tuple[1]
    linestyle = None
    color = None
    method_marker = "."
    markersize = 12
    plt.plot(param_name, method, data=df, marker=method_marker, linewidth=2, label = method_alias, markersize=markersize, linestyle = linestyle, color = color)
  plt.legend()
  if (xlabel != None):
    plt.xlabel(xlabel)
  if (ylabel != None):
    plt.ylabel(ylabel)
  plt.savefig(output)
  print("Saving result in " + output)
  plt.close()


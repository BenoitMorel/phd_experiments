import sys
import os
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
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

def get_methods_values(dataset, metric, methods, subst_model):
  values = {}
  print(dataset)
  for method in methods:
    values[method] = 0.0
    key = (method + "." + subst_model).lower()
    dataset_dir = os.path.join(exp.families_datasets_root, dataset)
    metrics = saved_metrics.get_metrics(dataset_dir, metric)
    print(metrics)
    values[method] = float(metrics[key])
  return values


def plot():
  output = "ensembl105_maxgapratio.svg"
  base = "ensembl105_single"
  max_gaps = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
  method_tuples = []
  method_tuples.append(("astrid-fastme_raxml-ng", "ASTRID"))
  method_tuples.append(("astral_raxml-ng", "ASTRAL-III"))
  method_tuples.append(("aster_raxml-ng", "ASTER"))
  method_tuples.append(("fastrfs-raxml-ng_single", "FastRFS"))
  method_tuples.append(("asteroid-raxml-ng", "Asteroid"))
  methods = []
  for mt in method_tuples:
    methods.append(mt[0])
  subst_model = "GTR+G"
  param_name = "max_gap_ratio"
  param_to_datasets = [] 
  for max_gap in max_gaps:
    dataset = base
    if (max_gap != 1.0):
      dataset = dataset + "_maxgapratio" + str(max_gap)
    param_to_datasets.append((1.0 - max_gap, dataset))

  values = {}    
  for xvalue,xdatasets in param_to_datasets:
    values[xvalue] = get_methods_values(xdatasets, "species_unrooted_rf", methods, subst_model)
  df = _get_df(values, methods, param_name)
  for method_tuple in method_tuples:
    method = method_tuple[0]
    method_alias = method_tuple[1]
    linestyle = None
    color = None
    method_marker = "."
    markersize = 12
    print(param_name)
    print(method)
    print(df)
    plt.plot(param_name, method, data=df, marker=method_marker, linewidth=2, label = method_alias, markersize=markersize, linestyle = linestyle, color = color)
  plt.legend()
  plt.xlabel("Gap ratio threshold") 
  plt.ylabel("RF distance")
  plt.savefig(output)
  print("Saving result in " + output)
  plt.close()





if (__name__ == "__main__"):
  plot()  

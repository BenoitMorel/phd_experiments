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
import common

def plot(datasets_rf_dico, x_param, fixed_params_dico, methods, x_label, output):
  print(x_param)
  datasets_to_plot = common.get_datasets_to_plot(datasets_rf_dico, fixed_params_dico, x_param)
  print(datasets_to_plot)
  datasets_to_plot.sort(key = lambda t: float(common.get_param_from_dataset_name(x_param, t)))
  df = pd.DataFrame()
  f, ax = plt.subplots(1)
 
  
  fake_df = {}
  fake_df[x_param] = []
  for method in methods:
    fake_df[method] = []
  for dataset in datasets_to_plot:
    print(float(common.get_param_from_dataset_name(x_param,     dataset)))
    fake_df[x_param].append(float(common.get_param_from_dataset_name(x_param, dataset)))
    rf_dico = datasets_rf_dico[dataset]
    for method in methods:
      if (not method in rf_dico):
        print("Warning: missing data for method " + method + " and dataset " + dataset)
    for method in rf_dico:
      if (not method in methods):
        print("Unknown method " + method)
        continue
      fake_df[method].append(float(rf_dico[method]))
  for elem in fake_df:
    print(elem)
    print(fake_df[elem])
    df[elem] = fake_df[elem]
  
  for method in methods:
    style = "solid"
    plt.plot(x_param, method, data=df, marker='x', linestyle = style, linewidth=2, label = method)
    ax.set_ylim(bottom=0)
  plt.xlabel(x_label[x_param])
  plt.ylabel('RF distance')
  plt.legend()
  plt.savefig(output)
  print("Saving result in " + output)
  plt.close()

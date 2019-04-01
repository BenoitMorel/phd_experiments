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

def plot(datasets_rf_dico, x_param, fixed_params_dico, methods, x_label, y_label, output):
  print(x_param)
  datasets_to_plot = common.get_datasets_to_plot(datasets_rf_dico, fixed_params_dico, x_param)
  #print(datasets_to_plot)
  datasets_to_plot.sort(key = lambda t: float(common.get_param_from_dataset_name(x_param, t)))
  df = pd.DataFrame()
  f, ax = plt.subplots(1)
 
  
  fake_df = {}
  fake_df[x_param] = []
  for method in methods:
    fake_df[method] = []
  for dataset in datasets_to_plot:
    #print(float(common.get_param_from_dataset_name(x_param,     dataset)))
    fake_df[x_param].append(float(common.get_param_from_dataset_name(x_param, dataset)))
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


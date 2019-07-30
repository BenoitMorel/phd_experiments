import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/generax')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/plotters')
import experiments as exp
import plot_histogram
import common
import scaling_generax
import fam
import boxplot
import rf_cells
import saved_metrics

import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style("darkgrid")

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


def plot(datasets_rf_dico, x_param, fixed_params_dico, methods, methods_dico, x_label, y_label, output):
  datasets_to_plot = get_datasets_to_plot(datasets_rf_dico, fixed_params_dico, x_param)
  datasets_to_plot.sort(key = lambda t: float(fam.get_param_from_dataset_name(x_param, t)))
  df = pd.DataFrame()
  f, ax = plt.subplots(1)
  fake_df = {}
  fake_df[x_param] = []
  for method in methods:
    fake_df[method] = []
  for dataset in datasets_to_plot:
    fake_df[x_param].append(float(fam.get_param_from_dataset_name(x_param, dataset)))
    rf_dico = datasets_rf_dico[dataset]
    for method in methods:
      method_pair = "true.true - " + method
      if (not method_pair in rf_dico):
        print("Warning: missing data for method " + method + " and dataset " + dataset)
    for method_pair in rf_dico:
      method = method_pair.split(" - ")[1]
      if (not method in methods):
        continue
      fake_df[method].append(float(rf_dico[method_pair]))
  for elem in fake_df:
    df[elem] = fake_df[elem]
 
  for method in methods:
    style = "solid"
    plt.plot(x_param, method, data=df, marker='.', linestyle = style, linewidth=2, label = methods_dico[method], markersize=12)
    #ax.set_ylim(bottom=0)
  plt.xlabel(x_label[x_param])
  plt.ylabel(y_label)
  plt.legend()
  plt.savefig(output)
  print("Saving result in " + output)
  plt.close()


class Plotter(object):
  def __init__(self, datasets_values_dico, fixed_parameters, methods, methods_dico, x_labels, y_label, prefix):
    self.datasets_values_dico = datasets_values_dico
    self.fixed_parameters = fixed_parameters
    self.methods = methods
    self.methods_dico = methods_dico
    self.x_labels = x_labels
    self.prefix = prefix
    self.y_label = y_label
  
  def __call__(self, parameter):
    plot(self.datasets_values_dico, parameter, self.fixed_parameters, self.methods, self.methods_dico, self.x_labels,  self.y_label, self.prefix + "_" + parameter + ".svg")



def get_runs(methods, model):
  return [fam.get_run_name(method, model) for method in methods]


def keys_sorted_by_values(dico, all_methods):
  s = sorted(dico, key=dico.get) 
  res = []
  for k in s:
    key = k.split(" - ")[1]
    if (key in all_methods and not "true.true" in key):
      res.append(key)
  return res
#return sorted(dico.items(), key=lambda kv: kv[1])

def compute_best_method_percentage(datasets_rf_dico_dl, datasets_rf_dico_dtl, subst_model, best_methods_to_quantify, all_methods):
  methods_set = set([])
  for method in best_methods_to_quantify:
    methods_set.add(method + "." + subst_model)
  total_it = 0
  best_it = 0
  dicos = [datasets_rf_dico_dl, datasets_rf_dico_dtl]
  for dico in dicos:
    for run in dico:
      dico_run = dico[run]
      sorted_methods = keys_sorted_by_values(dico_run, all_methods)
      if (sorted_methods[0] in methods_set):
        best_it += 1
      else:
        print(sorted_methods[0] + " is best")
      total_it += 1
  percentage = 100.0 * float(best_it) / float(total_it) 
  print("Methods " + str(best_methods_to_quantify) + " are best in " + str(percentage) + "% of the cases")

def plot_rrf(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl, subst_model):
  methods_dico = {}
  methods_dico["raxml-ng." + subst_model] = "RAxML-NG"
  methods_dico["notung90." + subst_model] = "Notung"
  methods_dico["phyldog." + subst_model] = "Phyldog"
  methods_dico["treerecs." + subst_model] = "Treerecs"
  methods_dico["ale-dl." + subst_model] = "ALE-DL"
  methods_dico["ale-dtl." + subst_model] = "ALE-DTL"
  methods_dico["generax-dl-random." + subst_model] = "GeneRax-DL"
  methods_dico["generax-dtl-random." + subst_model] = "GeneRax-DTL"
 
  

  methods_dl = get_runs(["raxml-ng", "notung90", "phyldog", "treerecs", "ale-dl", "generax-dl-random"], subst_model)
  methods_dtl = get_runs(["raxml-ng", "notung90", "phyldog", "treerecs", "ale-dtl", "generax-dtl-random"], subst_model)
  
  rf_y_label =  "Average relative RF"
  datasets_rf_dico_dl = common.get_metrics_for_datasets("jsim_", "average_rrf")
  datasets_rf_dico_dtl = common.get_metrics_for_datasets("jsimdtl_", "average_rrf")
  plotter_rf_dl = Plotter(datasets_rf_dico_dl, fixed_params_dl, methods_dl, methods_dico, x_labels, rf_y_label, "dl_rf")
  plotter_rf_dtl = Plotter(datasets_rf_dico_dtl, fixed_params_dtl, methods_dtl, methods_dico, x_labels, rf_y_label, "dtl_rf")
  for param in params_to_plot_dl:
    plotter_rf_dl(param)
  for param in params_to_plot_dtl:
    plotter_rf_dtl(param)
  
  #best_methods_to_quantify = ["ale-dl", "ale-dtl", "generax-dl-random", "generax-dtl-random"]
  best_methods_to_quantify = ["generax-dl-random", "generax-dtl-random"]
  compute_best_method_percentage(datasets_rf_dico_dl, datasets_rf_dico_dtl, subst_model, best_methods_to_quantify, methods_dico)

def plot_simulated_metrics():
  subst_model = "gtr+g"
    
  x_labels = {}
  x_labels["species"] = "Number of taxa in the species tree"
  x_labels["loss_rate"] = "Average D(T)L rates"
  x_labels["bl"] = "Gene tree branch length multiplier"
  x_labels["dl_ratio"] = "Rates ratio D/L (L is fixed)"
  x_labels["sites"] = "Number of sites"
  x_labels["dt_ratio"] = "Rates ratio D/T (D+T is fixed)"
  x_labels["perturbation"] = "Species tree relative RF distance to the true species tree"

  params_to_plot_dl = ["bl", "sites", "loss_rate", "dl_ratio", "species", "perturbation"]
  fixed_params_dl = {}
  fixed_params_dl["species"] = "19"
  fixed_params_dl["bl"] = "0.5"
  fixed_params_dl["loss_rate"] = "0.25"
  fixed_params_dl["dl_ratio"] = "1.0"
  fixed_params_dl["sites"] = "250"
  fixed_params_dl["families"] = "100"
  fixed_params_dl["perturbation"] = "0.0"
  
  params_to_plot_dtl = ["sites", "loss_rate", "bl", "species", "dt_ratio", "perturbation"]
  fixed_params_dtl = {}
  fixed_params_dtl["species"] = "19"
  fixed_params_dtl["loss_rate"] = "0.2"
  fixed_params_dtl["dt_ratio"] = "1.0"
  fixed_params_dtl["bl"] = "0.5"
  fixed_params_dtl["families"] = "100"
  fixed_params_dtl["sites"] = "250"
  fixed_params_dtl["perturbation"] = "0.0"
  
  plot_rrf(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl, subst_model)



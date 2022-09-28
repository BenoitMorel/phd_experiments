import os
import sys
import fam_data
import copy

# get a list of replicated datasets with varying parameters
# from a reference dataset, a list of parameters to vary, 
# and a range of replicates
def get_dataset_list_list(ref_dataset, strings_to_replace, replicates, add_ref = False):
  seeds_to_replace = []
  for i in replicates:
    seeds_to_replace.append("seed" + str(i))
  unreplicated_datasets = []
  if (add_ref):
    strings_to_replace = copy.deepcopy(strings_to_replace)
    strings_to_replace.append("none")
  fam_data.get_dataset_variations(unreplicated_datasets, ref_dataset, strings_to_replace)
  replicated_datasets = []
  for d in unreplicated_datasets:
    replicated = []
    fam_data.get_dataset_variations(replicated, d, seeds_to_replace)
    replicated_datasets.append(replicated)
  return replicated_datasets

def get_dataset_list(ref_dataset, strings_to_replace, replicates, add_ref = False):
  dataset_list_list = get_dataset_list_list(ref_dataset, strings_to_replace, replicates, add_ref)
  res = []
  for dataset_list in dataset_list_list:
    res.extend(dataset_list)
  return res

def get_plot_name(simulation_name, param, subst_model, metric_name):
  return "plot_" + simulation_name + "_" + metric_name.replace("_", "-") + "_" + param + "_" + subst_model + ".svg"



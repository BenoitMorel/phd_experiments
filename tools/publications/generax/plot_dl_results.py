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



def plot(datasets_rf_dico, x_param, fixed_params_dico, x_label, output):
  print(x_param)
  datasets_to_plot = common.get_datasets_to_plot(datasets_rf_dico, fixed_params_dico, x_param)
  print(datasets_to_plot)
  datasets_to_plot.sort(key = lambda t: float(common.get_param_from_dataset_name(x_param, t)))
  df = pd.DataFrame()
  f, ax = plt.subplots(1)
 
  
  methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Raxml", "GeneRax-DL-Random"]#, "GeneRax-DTL-Raxml", "GeneRax-DTL-Random"]
  fake_df = {}
  fake_df[x_param] = []
  for method in methods:
    fake_df[method] = []
  for dataset in datasets_to_plot:
    print(float(common.get_param_from_dataset_name(x_param,     dataset)))
    fake_df[x_param].append(float(common.get_param_from_dataset_name(x_param, dataset)))
    rf_dico = datasets_rf_dico[dataset]
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



if (__name__ == "__main__"): 

  datasets = common.get_available_datasets("jsim_")
  datasets_rf_dico = {}
  index = 0
  total = len(datasets)
  for dataset in datasets:
    print("Analyze " + dataset + " " + str(index) + "/" + str(total) )
    res = common.get_results(dataset)
    if (res != None):
      datasets_rf_dico[dataset] = res
    index += 1

  print("Plot...")

  x_label = {}
  x_label["species"] = "Number of taxa in the species tree"
  x_label["dup_rate"] = "Duplication rate"
  x_label["bl"] = "Gene tree branch length multiplier"
  x_label["dl_ratio"] = "Ratio between duplication and loss rates"
  x_label["sites"] = "Number of sites"

  
  default_species = "12"
  default_dup_rate = "0.5"
  default_bl = "1.0"
  default_dl_ratio = "2.0"
  default_sites = "500"
  default_families = "50"

  params_value_dico_sites = {}
  params_value_dico_sites["species"] = default_species
  params_value_dico_sites["dup_rate"] = default_dup_rate
  params_value_dico_sites["bl"] = default_bl
  params_value_dico_sites["dl_ratio"] = default_dl_ratio
  params_value_dico_sites["sites"] = default_sites
  params_value_dico_sites["families"] = default_families
  
  plot(datasets_rf_dico, "sites", params_value_dico_sites, x_label,  "sites.png")
  plot(datasets_rf_dico, "dup_rate", params_value_dico_sites, x_label, "rates.png")
  plot(datasets_rf_dico, "bl", params_value_dico_sites, x_label, "bl.png")
  plot(datasets_rf_dico, "species", params_value_dico_sites, x_label, "species.png")
  plot(datasets_rf_dico, "dl_ratio", params_value_dico_sites, x_label, "dl_ratio.png")

  default_families = "100"
  plot(datasets_rf_dico, "species", params_value_dico_sites, x_label, "species_2.png")
  




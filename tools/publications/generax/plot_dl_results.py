import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp
import common
import common_plots



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
  methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Raxml", "GeneRax-DL-Random"]#, "GeneRax-DTL-Raxml", "GeneRax-DTL-Random"]

  x_label = {}
  x_label["species"] = "Number of taxa in the species tree"
  x_label["dup_rate"] = "Duplication rate"
  x_label["bl"] = "Gene tree branch length multiplier"
  x_label["dl_ratio"] = "Ratio between duplication and loss rates"
  x_label["sites"] = "Number of sites"
  x_label["families"] = "Number of gene families"
  
  default_species = "19"
  default_dup_rate = "0.5"
  default_bl = "1.0"
  default_dl_ratio = "2.0"
  default_sites = "500"

  params_value_dico_sites = {}
  params_value_dico_sites["species"] = default_species
  params_value_dico_sites["dup_rate"] = default_dup_rate
  params_value_dico_sites["bl"] = default_bl
  params_value_dico_sites["dl_ratio"] = default_dl_ratio
  params_value_dico_sites["sites"] = default_sites
  params_value_dico_sites["families"] = "100"

  common_plots.plot(datasets_rf_dico, "sites", params_value_dico_sites, methods, x_label,  "sites.png")
  common_plots.plot(datasets_rf_dico, "dup_rate", params_value_dico_sites, methods, x_label, "rates.png")
  common_plots.plot(datasets_rf_dico, "bl", params_value_dico_sites, methods,  x_label, "bl.png")
  common_plots.plot(datasets_rf_dico, "species", params_value_dico_sites, methods, x_label, "species.png")
  common_plots.plot(datasets_rf_dico, "dl_ratio", params_value_dico_sites, methods, x_label, "dl_ratio.png")
  common_plots.plot(datasets_rf_dico, "families", params_value_dico_sites, methods, x_label, "families.png")

  

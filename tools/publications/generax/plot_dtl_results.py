import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp
import common
import common_plots


if (__name__ == "__main__"): 

  datasets = common.get_available_datasets("jsimdtl_")
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
  methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Random", "GeneRax-DTL-Random"]

  x_label = {}
  x_label["species"] = "Number of taxa in the species tree"
  x_label["dup_rate"] = "Average DTL rate"
  x_label["bl"] = "Gene tree branch length multiplier"
  x_label["dl_ratio"] = "Ratio between duplication and loss rates (loss rate is fixed)"
  x_label["sites"] = "Number of sites"
  x_label["tl_ratio"] = "Ratio between transfer and loss rates (loss rate is fixed"
  
  params_value_dico_sites = {}
  params_value_dico_sites["species"] = "16"
  params_value_dico_sites["dup_rate"] = "0.15"
  params_value_dico_sites["tl_ratio"] = "1.0"
  params_value_dico_sites["bl"] = "1.0"
  params_value_dico_sites["dl_ratio"] = "1.0"
  params_value_dico_sites["sites"] = "500"
  
  common_plots.plot(datasets_rf_dico, "sites", params_value_dico_sites, methods, x_label,  "sites_dtl.png")
  common_plots.plot(datasets_rf_dico, "species", params_value_dico_sites, methods, x_label,  "species_dtl.png")
  common_plots.plot(datasets_rf_dico, "tl_ratio", params_value_dico_sites, methods, x_label,  "transfers_dtl.png")
  common_plots.plot(datasets_rf_dico, "dup_rate", params_value_dico_sites, methods, x_label,  "dtl_rates_multiplier_dtl.png")







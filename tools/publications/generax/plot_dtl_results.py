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
  datasets_runtimes_dico = {}
  index = 0
  total = len(datasets)
  for dataset in datasets:
    print("Analyze " + dataset + " " + str(index) + "/" + str(total) )
    res = common.get_results(dataset)
    if (res != None):
      datasets_rf_dico[dataset] = res
    res = common.get_runtimes(dataset)
    if (res != None):
      datasets_runtimes_dico[dataset] = res
    index += 1

  print("Plot...")
  methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Random", "GeneRax-DTL-Random"]

  x_labels = {}
  x_labels["species"] = "Number of taxa in the species tree"
  x_labels["dup_rate"] = "Average DTL rate"
  x_labels["bl"] = "Gene tree branch length multiplier"
  x_labels["dl_ratio"] = "Ratio between duplication and loss rates (loss rate is fixed)"
  x_labels["sites"] = "Number of sites"
  x_labels["tl_ratio"] = "Ratio between transfer and loss rates (loss rate is fixed"
  
  params_value_dico_sites = {}
  params_value_dico_sites["species"] = "16"
  params_value_dico_sites["dup_rate"] = "0.15"
  params_value_dico_sites["tl_ratio"] = "1.0"
  params_value_dico_sites["bl"] = "1.0"
  params_value_dico_sites["dl_ratio"] = "1.0"
  params_value_dico_sites["sites"] = "500"

  plotter = common_plots.Plotter(datasets_rf_dico, datasets_runtimes_dico, params_value_dico_sites, methods, x_labels, "dtl")
  
  plotter("sites")
  plotter("species")
  plotter("tl_ratio")
  plotter("dup_rate")
  plotter("bl")







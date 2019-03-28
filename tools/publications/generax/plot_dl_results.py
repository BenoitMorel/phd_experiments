import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp
import common
import common_plots


if (__name__ == "__main__"): 


  x_labels = {}
  x_labels["species"] = "Number of taxa in the species tree"
  x_labels["dup_rate"] = "Duplication rate"
  x_labels["bl"] = "Gene tree branch length multiplier"
  x_labels["dl_ratio"] = "Ratio between duplication and loss rates"
  x_labels["sites"] = "Number of sites"
  x_labels["families"] = "Number of gene families"
  x_labels["tl_ratio"] = "Ratio between transfer and loss rates (loss rate is fixed"
  
  methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Raxml", "GeneRax-DL-Random"]#, "GeneRax-DTL-Raxml", "GeneRax-DTL-Random"]
  prefix = "jsim_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
  params_value_dico_sites = {}
  params_value_dico_sites["species"] = "19"
  params_value_dico_sites["dup_rate"] = "0.5"
  params_value_dico_sites["bl"] = "1.0"
  params_value_dico_sites["dl_ratio"] = "2.0"
  params_value_dico_sites["sites"] = "500"
  params_value_dico_sites["families"] = "100"
  plotter = common_plots.Plotter(datasets_rf_dico, datasets_runtimes_dico, params_value_dico_sites, methods, x_labels, "dl")
  plotter("sites")
  plotter("dup_rate")
  plotter("bl")
  plotter("dl_ratio")
  plotter("species")
  plotter("families")

  methods = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Random", "GeneRax-DTL-Random"]
  prefix = "jsimdtl_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
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






  

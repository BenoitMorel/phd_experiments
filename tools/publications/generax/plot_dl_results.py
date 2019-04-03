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
  
  #methods_dl_rf = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Raxml", "GeneRax-DL-Random"]
  #methods_dl_runtimes = ["RAxML-light", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Raxml", "GeneRax-DL-Random"]
  #methods_dtl_rf = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Random", "GeneRax-DTL-Random"]
  #methods_dtl_runtimes = ["RAxML-light", "Notung", "Phyldog", "Treerecs", "GeneRax-DL-Random", "GeneRax-DTL-Random"]
  
  methods_dl_rf = ["RAxML-NG", "ALE-DL", "GeneRax-DL-Raxml"]
  methods_dl_runtimes = list(methods_dl_rf)
  methods_dl_runtimes[0] = "RAxML-light"
  methods_dtl_rf = ["RAxML-NG", "ALE-DTL", "GeneRax-DTL-Raxml"]
  methods_dtl_runtimes = list(methods_dtl_rf)
  methods_dtl_runtimes[0] = "RAxML-light"
  
  params_to_plot_dl = ["sites", "dup_rate", "bl", "dl_ratio", "species", "families"]
  params_to_plot_dtl = ["sites", "dup_rate", "bl", "species", "tl_ratio"]
  
  rf_y_label =  "Average relative RF"
  runtime_y_label = "Runtime with 40 cores (s)"

  prefix = "jsim_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
  fixed_params = {}
  fixed_params["species"] = "19"
  fixed_params["dup_rate"] = "0.5"
  fixed_params["bl"] = "1.0"
  fixed_params["dl_ratio"] = "2.0"
  fixed_params["sites"] = "500"
  fixed_params["families"] = "100"
  plotter_rf = common_plots.Plotter(datasets_rf_dico, fixed_params, methods_dl_rf, x_labels, rf_y_label, "dl_rf")
  plotter_runtimes = common_plots.Plotter(datasets_runtimes_dico, fixed_params, methods_dl_runtimes, x_labels, runtime_y_label, "dl_runtimes")
  for param in params_to_plot_dl:
    plotter_rf(param)
    plotter_runtimes(param)

  prefix = "jsimdtl_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
  fixed_params = {}
  fixed_params["species"] = "16"
  fixed_params["dup_rate"] = "0.15"
  fixed_params["tl_ratio"] = "1.0"
  fixed_params["bl"] = "1.0"
  fixed_params["dl_ratio"] = "1.0"
  fixed_params["sites"] = "500"
  plotter_rf = common_plots.Plotter(datasets_rf_dico, fixed_params, methods_dtl_rf, x_labels, rf_y_label, "dtl_rf")
  plotter_runtimes = common_plots.Plotter(datasets_runtimes_dico, fixed_params, methods_dtl_runtimes, x_labels, runtime_y_label, "dtl_runtimes")
  for param in params_to_plot_dtl:
    plotter_rf(param)
    plotter_runtimes(param)






  

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
  x_labels["loss_rate"] = "Average D(T)L rates"
  x_labels["bl"] = "Gene tree branch length multiplier"
  x_labels["dl_ratio"] = "Ratio between duplication and loss rates"
  x_labels["sites"] = "Number of sites"
  x_labels["families"] = "Number of gene families"
  x_labels["dt_ratio"] = "Ratio between duplication and transfer (their sum is fixed)"
  
  methods_dl_rf = ["RAxML-NG", "Treerecs", "ALE-DL", "GeneRax-DL-Random"]
  methods_dtl_rf = ["RAxML-NG", "Treerecs", "ALE-DTL", "GeneRax-DTL-Random"]

  methods_dl_runtimes = list(methods_dl_rf)
  methods_dl_runtimes[0] = "RAxML-light"
  methods_dtl_runtimes = list(methods_dtl_rf)
  methods_dtl_runtimes[0] = "RAxML-light"
  
  params_to_plot_dl = ["sites", "loss_rate", "bl", "dl_ratio", "species", "families"]
  params_to_plot_dtl = ["sites", "loss_rate", "bl", "species", "dt_ratio", "families"]
  
  rf_y_label =  "Average relative RF"
  runtime_y_label = "Runtime with 40 cores (s)"

  prefix = "jsim_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
  fixed_params = {}
  fixed_params["species"] = "19"
  fixed_params["bl"] = "0.5"
  fixed_params["loss_rate"] = "0.25"
  fixed_params["dl_ratio"] = "1.0"
  fixed_params["sites"] = "500"
  fixed_params["families"] = "100"
  
  plotter_rf = common_plots.Plotter(datasets_rf_dico, fixed_params, methods_dl_rf, x_labels, rf_y_label, "dl_rf")
  plotter_runtimes = common_plots.Plotter(datasets_runtimes_dico, fixed_params, methods_dl_runtimes, x_labels, runtime_y_label, "dl_runtimes")
  for param in params_to_plot_dl:
    plotter_rf(param)
    plotter_runtimes(param)

  
  #fixed_point_dtl = "jsimdtl_s16_f100_sites500_dna4_bl0.5_d0.1_l0.2_t0.1"

  prefix = "jsimdtl_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
  fixed_params = {}
  fixed_params["species"] = "16"
  fixed_params["loss_rate"] = "0.2"
  fixed_params["dt_ratio"] = "1.0"
  fixed_params["bl"] = "0.5"
  fixed_params["families"] = "100"
  fixed_params["sites"] = "500"
  plotter_rf = common_plots.Plotter(datasets_rf_dico, fixed_params, methods_dtl_rf, x_labels, rf_y_label, "dtl_rf")
  plotter_runtimes = common_plots.Plotter(datasets_runtimes_dico, fixed_params, methods_dtl_runtimes, x_labels, runtime_y_label, "dtl_runtimes")
  for param in params_to_plot_dtl:
    plotter_rf(param)
    plotter_runtimes(param)






  

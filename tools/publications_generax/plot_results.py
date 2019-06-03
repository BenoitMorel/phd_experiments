import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp
import common


def substract(datasets_dico, method):
  for dataset in datasets_dico:
    print("substract dataset " + dataset)
    dico = datasets_dico[dataset]
    substract = float(dico[method])
    for m in dico:
      dico[m] = str(float(dico[m]) - substract)
    print(dataset)
    print(datasets_dico[dataset])

def plot_rrf(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl):
  methods_dl = ["raxml-ng", "notung80", "phyldog", "treerecs", "ale-dl", "generax-dl-raxml"]
  methods_dtl = ["raxml-ng", "notung80", "phyldog", "treerecs", "ale-dtl", "generax-dtl-raxml"] 
  rf_y_label =  "Average relative RF"
  datasets_rf_dico_dl = common.get_metrics_for_datasets("jsim_", "average_rrf")
  datasets_rf_dico_dtl = common.get_metrics_for_datasets("jsimdtl_", "average_rrf")
  plotter_rf_dl = common.Plotter(datasets_rf_dico_dl, fixed_params_dl, methods_dl, x_labels, rf_y_label, "dl_rf")
  plotter_rf_dtl = common.Plotter(datasets_rf_dico_dtl, fixed_params_dtl, methods_dl, x_labels, rf_y_label, "dtl_rf")
  for param in params_to_plot_dl:
    plotter_rf_dl(param)
  for param in params_to_plot_dtl:
    plotter_rf_dtl(param)


def plot_runtimes(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl):
  methods_dl = ["raxml-light", "notung80", "treerecs", "ale-dl", "generax-dl-raxml"] 
  methods_dtl = ["raxml-light", "notung80", "treerecs", "ale-dtl", "generax-dtl-raxml"]
  runtime_y_label = "Runtime with 40 cores (s)"
  datasets_dico_dl = common.get_metrics_for_datasets("jsim_", "runtimes")
  datasets_dico_dtl = common.get_metrics_for_datasets("jsimdtl_", "runtimes")
  plotter_runtimes_dl = common.Plotter(datasets_dico_dl, fixed_params_dl, methods_dl, x_labels, runtime_y_label, "dl_runtimes")
  plotter_runtimes_dtl = common.Plotter(datasets_dico_dtl, fixed_params_dtl, methods_dtl, x_labels, runtime_y_label, "dtl_runtimes")
  for param in params_to_plot_dl:
    plotter_runtimes_dl(param)
  for param in params_to_plot_dtl:
    plotter_runtimes_dtl(param)

def plot_ll(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl):
  methods_dl_rf = ["true", "raxml-ng", "notung80", "phyldog", "treerecs", "ale-dl", "generax-dl-raxml"]
  methods_dtl_rf = ["true", "raxml-ng", "notung80", "phyldog", "treerecs", "ale-dtl", "generax-dtl-raxml"] 
  ll_DL_y_label = "Joint likelihood difference with true trees"
  ll_DTL_y_label = "Joint likelihood difference with true trees"
  datasets_dico_dl = common.get_metrics_for_datasets("jsim_", "JointLL_DL")
  datasets_dico_dtl = common.get_metrics_for_datasets("jsimdtl_", "JointLL_DTL")
  try:
    substract(datasets_dico_dl, "True")
    substract(datasets_dico_dtl, "True")
  except:
    print("failed substracting likelihoods")
    return
  plotter_dl = common.Plotter(datasets_dico_dl, fixed_params_dl, methods_dl_ll, x_labels, ll_DL_y_label, "dl_ll")
  plotter_dtl = common.Plotter(datasets_dico_dtl, fixed_params_dtl, methods_dtl_ll, x_labels, ll_DTL_y_label, "dtl_ll")
  for param in params_to_plot_dl:
    plotter_runtimes_dl(param)
  for param in params_to_plot_dtl:
    plotter_runtimes_dtl(param)

if (__name__ == "__main__"): 
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
  fixed_params_dl["sites"] = "500"
  fixed_params_dl["families"] = "100"
  fixed_params_dl["perturbation"] = "0.0"
  
  params_to_plot_dtl = ["sites", "loss_rate", "bl", "species", "dt_ratio", "perturbation"]
  fixed_params_dtl = {}
  fixed_params_dtl["species"] = "19"
  fixed_params_dtl["loss_rate"] = "0.2"
  fixed_params_dtl["dt_ratio"] = "1.0"
  fixed_params_dtl["bl"] = "0.5"
  fixed_params_dtl["families"] = "100"
  fixed_params_dtl["sites"] = "500"
  fixed_params_dtl["perturbation"] = "0.0"
  
  plot_runtimes(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl)
  plot_rrf(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl)
  plot_ll(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl)



  

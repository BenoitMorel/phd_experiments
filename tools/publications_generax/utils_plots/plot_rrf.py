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



def get_runs(methods, model):
  return [fam.get_run_name(method, model) for method in methods]

def plot_rrf(x_labels, params_to_plot_dl, params_to_plot_dtl, fixed_params_dl, fixed_params_dtl, subst_model):
  methods_dl = get_runs(["raxml-ng", "notung80", "phyldog", "treerecs", "ale-dl", "generax-dl-raxml"], subst_model)
  methods_dtl = get_runs(["raxml-ng", "notung80", "phyldog", "treerecs", "ale-dtl", "generax-dtl-raxml"], subst_model)
  rf_y_label =  "Average relative RF"
  datasets_rf_dico_dl = common.get_metrics_for_datasets("jsim_", "average_rrf")
  datasets_rf_dico_dtl = common.get_metrics_for_datasets("jsimdtl_", "average_rrf")
  plotter_rf_dl = common.Plotter(datasets_rf_dico_dl, fixed_params_dl, methods_dl, x_labels, rf_y_label, "dl_rf")
  plotter_rf_dtl = common.Plotter(datasets_rf_dico_dtl, fixed_params_dtl, methods_dtl, x_labels, rf_y_label, "dtl_rf")
  for param in params_to_plot_dl:
    plotter_rf_dl(param)
  for param in params_to_plot_dtl:
    plotter_rf_dtl(param)

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


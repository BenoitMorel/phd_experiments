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

if (__name__ == "__main__"): 


  x_labels = {}
  x_labels["species"] = "Number of taxa in the species tree"
  x_labels["loss_rate"] = "Average D(T)L rates"
  x_labels["bl"] = "Gene tree branch length multiplier"
  x_labels["dl_ratio"] = "Rates ratio D/L (L is fixed)"
  x_labels["sites"] = "Number of sites"
  #x_labels["families"] = "Number of gene families"
  x_labels["dt_ratio"] = "Rates ratio D/T (D+T is fixed)"
  
  methods_dl_rf = ["RAxML-NG", "Notung", "Phyldog", "Treerecs", "ALE-DL", "GeneRax-DL-Random", "GeneRax-DL-Raxml"]
  methods_dtl_rf = ["RAxML-NG", "Treerecs", "ALE-DTL", "GeneRax-DTL-Random", "GeneRax-DTL-Raxml"]

  methods_dl_runtimes = list(methods_dl_rf)
  methods_dl_runtimes[0] = "RAxML-light"
  methods_dtl_runtimes = list(methods_dtl_rf)
  methods_dtl_runtimes[0] = "RAxML-light"
  
  methods_dl_ll = ["True"]
  methods_dl_ll.extend(methods_dl_rf)
  methods_dtl_ll = ["True"]
  methods_dtl_ll.extend(methods_dtl_rf)
  
  params_to_plot_dl = ["bl"] #, "sites", "loss_rate", "dl_ratio", "species"]
  params_to_plot_dtl = [] #"sites", "loss_rate", "bl", "species", "dt_ratio"]
  
  rf_y_label =  "Average relative RF"
  runtime_y_label = "Runtime with 40 cores (s)"
  ll_DL_y_label = "Joint likelihood difference with true trees"
  ll_DTL_y_label = "Joint likelihood difference with true trees"

  prefix = "jsim_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
  datasets_ll_dico = common.get_metrics_for_datasets(prefix, "joint_likelihood_DL")
  try:
    substract(datasets_ll_dico, "True")
  except:
    print("failed substracting likelihoods")

  fixed_params = {}
  fixed_params["species"] = "19"
  fixed_params["bl"] = "0.5"
  fixed_params["loss_rate"] = "0.25"
  fixed_params["dl_ratio"] = "1.0"
  fixed_params["sites"] = "500"
  fixed_params["families"] = "100"
  
  plotter_rf = common.Plotter(datasets_rf_dico, fixed_params, methods_dl_rf, x_labels, rf_y_label, "dl_rf")
  plotter_ll = common.Plotter(datasets_ll_dico, fixed_params, methods_dl_ll, x_labels, ll_DL_y_label, "dl_ll")
  plotter_runtimes = common.Plotter(datasets_runtimes_dico, fixed_params, methods_dl_runtimes, x_labels, runtime_y_label, "dl_runtimes")
  for param in params_to_plot_dl:
    plotter_rf(param)
    plotter_runtimes(param)
    #plotter_ll(param)
  
  #fixed_point_dtl = "jsimdtl_s16_f100_sites500_dna4_bl0.5_d0.1_l0.2_t0.1"

  prefix = "jsimdtl_"
  datasets_rf_dico = common.get_metrics_for_datasets(prefix, "average_rrf")
  datasets_runtimes_dico = common.get_metrics_for_datasets(prefix, "runtimes")
  datasets_ll_dico = common.get_metrics_for_datasets(prefix, "joint_likelihood_DTL")
  substract(datasets_ll_dico, "True")
  
  fixed_params = {}
  fixed_params["species"] = "16"
  fixed_params["loss_rate"] = "0.2"
  fixed_params["dt_ratio"] = "1.0"
  fixed_params["bl"] = "0.5"
  fixed_params["families"] = "100"
  fixed_params["sites"] = "500"
  plotter_rf = common.Plotter(datasets_rf_dico, fixed_params, methods_dtl_rf, x_labels, rf_y_label, "dtl_rf")
  plotter_ll = common.Plotter(datasets_ll_dico, fixed_params, methods_dtl_ll, x_labels, ll_DTL_y_label, "dtl_ll")
  plotter_runtimes = common.Plotter(datasets_runtimes_dico, fixed_params, methods_dtl_runtimes, x_labels, runtime_y_label, "dtl_runtimes")
  for param in params_to_plot_dtl:
    plotter_rf(param)
    plotter_runtimes(param)
    plotter_ll(param)





  

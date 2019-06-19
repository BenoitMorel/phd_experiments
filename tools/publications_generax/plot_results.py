import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/generax')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/plotters')
import experiments as exp
import common
import scaling_generax
import fam
import boxplot
import rf_cells

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

def plot_simulated_metrics():
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

def plot_scaling():
  datasets = ["ensembl_96_ncrna_primates", "cyano_empirical"]
  for dataset in datasets:
    datadir = fam.get_datadir(dataset)
    scaling_generax.plot_scaling_metric(datadir)

def plot_boxplots():
  methods = ["raxml", "treerecs", "notung80", "phyldog", "ale-dtl"]
  models = ["GTR+G"]
  datasets = ["cyano_simulated"]
  for model in models:
    print("model " + model)
    runs = []
    for method in methods:
      runs.append(fam.get_run_name(method, model))
    bp = boxplot.BoxPlot(ylabel = "rf distances")
    data = {}
    for run in runs:
      data[run] = []
    for dataset in datasets:
      cells = rf_cells.load_rf_cells(fam.get_datadir(dataset))
      for family in cells:
        family_cells = cells[family]
        for run in runs:
          cell = rf_cells.get_rf_to_true(family_cells, run)
          data[run].append(cell[0] / cell[1])
    for run in data:
      bp.add_elem(run, data[run])
    output = "plop." + model + ".svg"
    bp.plot(output)
    print("plot " + output)

def plot_model_boxplots():
  methods = ["raxml-ng", "treerecs", "notung80", "generax-dl-random"]
  models = ["JC", "GTR+G"]
  dico = {}
  datasets = ["jsimdtl_s5_f10_sites100_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0"]
  for model in models:
    model_dico = {}
    for method in methods:
      model_dico[method] = []
    dico[model] = model_dico

  for dataset in datasets:
    cells = rf_cells.load_rf_cells(fam.get_datadir(dataset))
    for family in cells:
      family_cells = cells[family]
      for model in models:
        for method in methods:
          cell = rf_cells.get_rf_to_true(family_cells, fam.get_run_name(method, model))
          dico[model][method].append(cell[0] / cell[1])
  gbp = boxplot.GroupBoxPlot(data = dico, title = None, xlabel = "Methods", ylabel = "Relative RF distance", hue_label = "Substitution model")
  output = os.path.abspath("grouped_box_plot.svg")
  gbp.plot(output)
  print("Output in " + output)



if (__name__ == "__main__"):
  #plot_simulated_metrics()
  #plot_scaling()  
  #plot_boxplots()
  plot_model_boxplots()


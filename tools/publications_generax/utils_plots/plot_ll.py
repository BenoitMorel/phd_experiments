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

def aux(dataset, model, recModel):
  methods_to_plot = []
  methods_to_plot.append("raxml-ng")
  methods_to_plot.append("notung90")
  methods_to_plot.append("treerecs")
  methods_to_plot.append("phyldog")
  methods_to_plot.append("eccetera")
  methods_to_plot.append("ale-" + recModel.lower())
  methods_to_plot.append("generax-" + recModel.lower() + "-random")

  methods_display_name = {}
  methods_display_name["raxml-ng"] = "RAxML-NG"
  methods_display_name["ale-" + recModel.lower()] = "ALE-" + recModel
  methods_display_name["generax-" + recModel.lower() + "-random"] = "GeneRax-" + recModel
  methods_display_name["treerecs"] = "Treerecs"
  methods_display_name["notung90"] = "Notung"
  methods_display_name["phyldog"] = "Phyldog"
  methods_display_name["eccetera"] = "EcceTERA"

  categories = {}
  categories["JointLL_DL"] = "DL"
  categories["JointLL_DTL"] = "DTL"

  datadir = fam.get_datadir(dataset)
  dico = {} # dico[cat][run] = runtime
  for cat in categories:
    dico[cat] = {}
  runs = saved_metrics.get_metrics_methods(datadir, "JointLL_DL")
  saved_metrics_dict = {}
  for cat in categories:
    saved_metrics_dict[cat] = saved_metrics.get_metrics(datadir, cat)
   
  for cat in categories:
    for run in runs:
      if (model.lower() in run.lower()):
        method = fam.get_method_from_run(run)
        dico[cat][method] = float(saved_metrics_dict[cat][run])
  
  x_order = []
  yvalues = {}
  for method in methods_to_plot:
    x_order.append(methods_display_name[method])
  for cat in categories:
    yvalues[categories[cat]] = {}
  for cat in categories:
    for method in methods_to_plot:
      x = methods_display_name[method]
      yvalues[categories[cat]][x] = dico[cat][method]
  
  output = "joint_likelihood__" + dataset
  # GROUPED HISTOGRAM
  plot_histogram.plot_grouped_histogram(yvalues, cat_name = "Category", class_name = "Methods", values_name = "Joint log-likelihood", start_at_min_y = True, x_order = x_order, output = output + ".svg")
  # SIMPLE HISTOGRAMS
  simple_xlabels = x_order
  simple_yvalues  = []
  for x in x_order:
    simple_yvalues.append(yvalues[recModel][x])
  output = "simple_ll_" + dataset + "_" + recModel + ".svg"
  plot_histogram.plot_histogram(simple_xlabels, simple_yvalues, ycaption = "Joint log-likelihood", start_at_min_y = True, reverse_bar = True,  output = output)

def plot_ll():
  datasets = ["cyano_empirical", "ensembl_96_ncrna_primates"]
  aux("cyano_empirical", "LG+G", "DTL")
  aux("ensembl_96_ncrna_primates", "GTR+G", "DL")

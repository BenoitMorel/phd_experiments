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


def plot(yvalues, output):
  print(yvalues)
  plot_histogram.plot_grouped_histogram(yvalues, cat_name = "Category", class_name = "Methods", values_name = "Time (s)", log_scale = True,  output = output)


def plot_runtimes():
  datasets = ["cyano_empirical"]
  model = "LG+G"
  methods_to_plot = []
  methods_to_plot.append(["raxml-light"])
  methods_to_plot.append(["notung80", "raxml-ng"])
  methods_to_plot.append(["treerecs", "raxml-ng"])
  #methods_to_plot.append(["phyldog", "raxml-ng"])
  methods_to_plot.append(["ale-dtl", "mrbayes"])
  methods_to_plot.append(["generax-dtl-raxml", "raxml-light"])
  methods_to_plot.append(["generax-dtl-random"])

  methods_display_name = {}
  methods_display_name["raxml-light"] = "RAxML-NG"
  methods_display_name["ale-dtl"] = "ALE"
  methods_display_name["generax-dtl-raxml"] = "GeneRax-raxml"
  methods_display_name["generax-dtl-random"] = "GeneRax-rdm"
  methods_display_name["treerecs"] = "Treerecs"
  methods_display_name["notung80"] = "Notung"

  for dataset in datasets:
    output = "runtimes_" + dataset + ".svg"
    categories = ["Parallel", "Sequential"]
    datadir = fam.get_datadir(dataset)
    dico = {} # dico[cat][run] = runtime
    for cat in categories:
      dico[cat] = {}
    runs = saved_metrics.get_metrics_methods(datadir, "runtimes")
    saved_metrics_runtimes = saved_metrics.get_metrics(datadir, "runtimes")
    saved_metrics_seqtimes = saved_metrics.get_metrics(datadir, "seqtimes")
    for run in runs:
      if (model.lower() in run.lower()):
        method = fam.get_method_from_run(run)
        dico["Parallel"][method] = float(saved_metrics_runtimes[run])
        if (run in saved_metrics_seqtimes):
          dico["Sequential"][method] = float(saved_metrics_seqtimes[run])
        else:
          dico["Sequential"][method] = float(saved_metrics_runtimes[run])
    xlabels = []
    yvalues_individual = {}
    yvalues_cumulated = {}
    for cat in categories:
      yvalues_individual[cat] = {}
      yvalues_cumulated[cat] = {}
    for cat in categories:
      for method in methods_to_plot:
        x = methods_display_name[method[0]]
        yvalues_individual[cat][x] = dico[cat][method[0]]
        cumulated_time = 0.0
        for m in method:
          cumulated_time += dico[cat][m]
        yvalues_cumulated[cat][x] = cumulated_time
    
    output = "runtimes_" + dataset
    plot(yvalues_individual, output + "_individual.svg")
    plot(yvalues_cumulated, output + "_cumulated.svg")
     


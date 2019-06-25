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
import math



def plot_boxplots():
  methods = ["raxml", "treerecs", "notung80", "phyldog", "ale-dtl"]
  models = ["gtr+g"]
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
  methods = ["raxml-ng", "ale-dtl", "generax-dtl-random"]
  models = ["LG+G", "DAYHOFF"]
  dico = {}
  datasets = ["cyano_simulated"]
  title = "Simulated cyanobacteriai dataset with 1099 families"
        
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
  gbp = boxplot.GroupBoxPlot(data = dico, title = title, ylabel = "Relative RF distance", hue_label = "Substitution model", order = methods)
  output = os.path.abspath("grouped_box_plot.svg")
  gbp.plot(output)
  print("Output in " + output)

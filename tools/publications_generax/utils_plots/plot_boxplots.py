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



def plot_model_boxplots():
  methods = {}
  methods["generax-dtl-random"] = "GeneRax"
  methods["ale-dtl"] = "ALE"
  #methods["deleterious"] = "JPrIME-DLTRS"
  #methods["treerecs"] = "Treerecs"
  methods["notung90"] = "Notung"
  methods["eccetera"] = "EcceTERA"
  #methods["phyldog"] = "Phyldog"
  methods["raxml-ng"] = "RAxML-NG"
  #order = ["GeneRax", "ALE", "JPrIME-DLTRS", "Treerecs",  "EcceTERA", "Phyldog", "Notung", "RAxML-NG"]
  order = ["GeneRax", "ALE", "EcceTERA", "Notung", "RAxML-NG"]
  models = ["LG+G+I", "WAG"]
  models_to_display = {}
  models_to_display["LG+G+I"] = "True model"
  models_to_display["WAG"] = "Wrong model"
  
  dico = {}
  datasets = ["cyano_simulated"]
  title = "Simulated cyanobacteria (1099 families)"
        
  for model in models:
    model_dico = {}
    for method in methods:
      model_dico[methods[method]] = []
    dico[models_to_display[model]] = model_dico

  for dataset in datasets:
    cells = rf_cells.load_rf_cells(fam.get_datadir(dataset))
    for family in cells:
      family_cells = cells[family]
      for model in models:
        for method in methods:
          cell = rf_cells.get_rf_to_true(family_cells, fam.get_run_name(method, model))
          dico[models_to_display[model]][methods[method]].append(cell[0] / cell[1])
  gbp = boxplot.GroupBoxPlot(data = dico, title = title, ylabel = "Relative RF distance", hue_label = "Substitution model", order = order)
  output = os.path.abspath("cyano_simulated_boxplot.svg")
  gbp.plot(output)
  print("Output in " + output)


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
import kf_cells
import saved_metrics
import math

def aux_rf(methods, order, models, models_to_display, datasets, title, output):
  dico = {}
        
  for model in models:
    model_dico = {}
    for method in methods:
      if (methods[method] in order):
        model_dico[methods[method]] = []
    dico[models_to_display[model]] = model_dico

  for dataset in datasets:
    cells = rf_cells.load_rf_cells(fam.get_datadir(dataset))
    for family in cells:
      family_cells = cells[family]
      for model in models:
        for method in methods:
          if (methods[method] in order):
            cell = rf_cells.get_rf_to_true(family_cells, fam.get_run_name(method, model))
            dico[models_to_display[model]][methods[method]].append(cell[0] / cell[1])
  gbp = boxplot.GroupBoxPlot(data = dico, title = title, ylabel = "Relative RF distance", hue_label = "Substitution model", order = order)
  gbp.plot(output)
  print("Output in " + output)



def aux_kf(methods, order, model, models_to_display, datasets, title, output, max_y):
        
  model_dico = {}
  for method in methods:
    if (methods[method] in order):
      model_dico[methods[method]] = []

  for dataset in datasets:
    cells = kf_cells.load_kf_cells(fam.get_datadir(dataset))

    for family in cells:
      family_cells = cells[family]
      for method in methods:
        if (methods[method] in order):
          cell = kf_cells.get_kf_to_true(family_cells, fam.get_run_name(method, model))
          model_dico[methods[method]].append(cell[0] / cell[1])
  bp = boxplot.BoxPlot(title, ylabel = "Branch Score Distance", order = order, max_y = max_y)
  for method in methods:
    if (methods[method] in order):
      bp.add_elem(methods[method], model_dico[methods[method]])
  bp.plot(output)
  print("Output in " + output)



def plot_model_boxplots():
  methods = {}
  #methods["generax_spr10-dtl-random"] = "GeneRax-spr10"
  methods["generax-dtl-random"] = "GeneRax"
  methods["ale-dtl"] = "ALE"
  #methods["deleterious"] = "JPrIME-DLTRS"
  methods["treerecs"] = "Treerecs"
  methods["notung90"] = "Notung"
  methods["eccetera"] = "EcceTERA"
  methods["phyldog"] = "Phyldog"
  methods["raxml-ng"] = "RAxML-NG"
  order = ["GeneRax", "ALE", "EcceTERA", "Treerecs", "Phyldog", "RAxML-NG", "Notung"]
  models = ["LG+G+I", "WAG"]
  models_to_display = {}
  models_to_display["LG+G+I"] = "True model"
  models_to_display["WAG"] = "Wrong model"
  datasets = ["cyano_simulated"]
  title = "Simulated cyanobacteria (1099 families)"
  output = os.path.abspath("cyano_simulated_boxplot.svg")
  aux_rf(methods, order, models, models_to_display, datasets, title, output)
  output = os.path.abspath("cyano_simulated_boxplot_kf.svg")
  model = "LG+G+I"
  order = ["GeneRax", "RAxML-NG", "Treerecs", "ALE", "Phyldog"]
# scale to 12, otherwise the scaling is to big because of phyldog
  aux_kf(methods, order, model, models_to_display, datasets, title, output, 12)
  
  





  

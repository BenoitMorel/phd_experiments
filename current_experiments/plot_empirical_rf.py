import os
import sys
sys.path.insert(0, os.path.join("tools", "families"))
import saved_metrics
import fam
sys.path.insert(0, os.path.join("tools", "plotters"))
import plot_histogram

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set_style("darkgrid")





# tuples of dataset_name and datatype 
# (0 for DNA, 1 for Amino Acids)
datasets = []
datasets.append(("cyano_empirical", 1))
datasets.append(("aa_ensembl_98_ncrna_primates", 0))
datasets.append(("aa_ensembl_98_ncrna_lowprimates", 0))
datasets.append(("aa_ensembl_98_ncrna_mammals", 0))
datasets.append(("aa_ensembl_98_ncrna_vertebrates", 0))
datasets.append(("aa_ensembl_98_ncrna_allvertebrates", 0))
metric_names = ["species_unrooted_rf"]

# methods_dict[formatted_name] = name_to_display
methods_dict = {}
methods_dict["astralpro"] = "AstralPro"
methods_dict["fastmulrfs-single"] = "FastMULRFS"
methods_dict["generax-MiniNJ-prune-fam"] = "GeneRaxMiniNJ"
methods_dict["duptree"] = "DupTree"
methods_dict["njrax-Cherry"] = "CherryMerging"
methods_dict["njrax-MiniNJ"] = "MiniNJ"

# dict of tuples: gene_trees[gene_name] = (dna_model, aa_model)
gene_trees = {}
gene_trees["true"] = ("true", "true")
gene_trees["raxml-ng"] = ("GTR+G", "LG+G+I")
gene_trees["fasttree"] = ("GTR", "LG")

def plot_grouped_histogram(data, title = None, xcaption = None, ycaption = None, x_order = None, output = "show"):
  values_vec = []
  categories_vec = []
  classes_vec = []
  min_value = 99999999999999
  max_value = -9999999999999
  cat_name = "categories"
  class_name = "classes"
  values_name = "values"
  kind = "bar"
  for category in data:
    for classs in data[category]:
      value = data[category][classs]
      min_value = min(min_value, value)
      max_value = max(max_value, value)
      values_vec.append(value)
      categories_vec.append(category)
      classes_vec.append(classs)
  plt.xticks(rotation = 45)
  dataFrameDico = {}
  dataFrameDico[class_name] = classes_vec
  dataFrameDico[cat_name] = categories_vec
  dataFrameDico[values_name] = values_vec
  df = pd.DataFrame(data = dataFrameDico)
  f = sns.factorplot(data = df, x = class_name, hue = cat_name, y = values_name, kind = kind, order = x_order, legend = False)
  if (kind == "strip"):
    plt.setp(f.ax.collections, sizes=[100])
  f.despine(left=True)
  #f._legend.set_title('')
  plt.gca().legend().set_title('')
  plt.gca().set_xlabel('')
  f.set(ylim=(-0.01, max_value * 2.0))
  plt.xticks(rotation = 45)
  if (title != None):
    f.fig.suptitle(title)
  f.fig.tight_layout()
  if (output == "show"):
    plt.show()
  else:
    plt.savefig(output)
  
    

    
def empirical_plot(dataset_name, dataset_type, metric_name, methods_dict, gene_trees):
  dataset_dir = fam.get_datadir(dataset_name)
  figure_name = "empirical_" + dataset_name + "_" + metric_name + ".svg"
  #data[category][class] = value
  # categories are gene_trees, and classes are species_methods 
  data = {} 
  # get all metric values for this dataset and metric_name
  metrics = saved_metrics.get_metrics(dataset_dir, metric_name)
  for gene_method in gene_trees:
    data[gene_method] = {}
    subst_model = gene_trees[gene_method][dataset_type] 
    for species_method in methods_dict:
      # get the value for this gene_method and species_method
      run_name = species_method + "_" + gene_method + "." + subst_model
      value = -1.0
      try:
        value = float(metrics[run_name.lower()])
      except:
        pass
      # store the value in the data dict
      data[gene_method][methods_dict[species_method]] = value
  print(data)
  title = dataset_name
  xcaption = "Species inference method"
  ycaption = "Relative RF distance"
  plot_grouped_histogram(data, title, xcaption, ycaption,  x_order = None, output = figure_name)
  print("Saved plot " + figure_name)

for dataset in datasets:
  dataset_name = dataset[0]
  dataset_type = dataset[1]
  for metric_name in metric_names:
    empirical_plot(dataset_name, dataset_type, metric_name, methods_dict, gene_trees)




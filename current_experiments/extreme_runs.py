import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "plotters"))
import fam
import fam_data
import saved_metrics
import run_all_species
import boxplot
from run_all_species import SpeciesRunFilter

datasets = []
cores = 40
launch_mode = "normal"

do_run = False
do_plot = True

run_replicates = range(2000, 2050)
plot_replicates = range(2000, 2050)
extreme_datasets = {}
sites_range = ["100", "200", "300"]
species_range = ["20", "25", "30"]
for sites in sites_range:
  for species in species_range:
# DTL
    extreme_datasets["DTL_s" + species + "_sites" + sites] = "ssim_var_s" + species + "_f100_sites" + sites + "_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop10_mu1.0_theta0.0_seed1000"
# ILS and DL
    extreme_datasets["ILS_DL_s" + species + "_sites" + sites] = "ssim_var_s" + species + "_f100_sites" + sites + "_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop100000000_mu1.0_theta0.0_seed1000"
# ILS and DTL
    #extreme_datasets["ILS_DTL_" + sites] = "ssim_var_s" + species + "_f100_sites" + sites + "_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop100000000_mu1.0_theta0.0_seed1000"



def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def get_dataset_list(ref_dataset, replicates):
  seeds_to_replace = []
  for i in replicates:
    seeds_to_replace.append("seed" + str(i))
  replicated_datasets = []
  fam_data.get_dataset_variations(replicated_datasets, ref_dataset, seeds_to_replace)
  return replicated_datasets

def run_extreme_cases():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.generate = True
  run_filter.pargenes_starting_trees = 1
  run_filter.pargenes_bootstrap_trees = 0
  run_filter.pargenes = True
  run_filter.duptree = True
  run_filter.njrax = True
  run_filter.astralpro = True
  run_filter.njst = True
  run_filter.cherry = True
  run_filter.speciesraxperfamily = True
  run_filter.cleanup = True
  #run_filter.astralpromultiple = True
  #run_filter.concatenation_min = True
  #run_filter.concatenation_max = True
  #run_filter.speciesrax = True
  datasets = []
  for dataset_name in extreme_datasets:
    dataset = extreme_datasets[dataset_name]
    datasets += get_dataset_list(dataset, run_replicates)
  run_species_methods(datasets, "GTR+G", cores, run_filter, launch_mode)


def get_method_to_values(datasets, methods, metric):
  res = {}
  for m in methods:
    res[m] = []
  for dataset in datasets:
    datadir = fam.get_datadir(dataset)
    metrics = saved_metrics.get_metrics(datadir, metric)
    for m in methods:
      if (not m in metrics):
        print("Cannot find metric " + metric + " for method " + m + " in " + dataset)
        return None
      res[m].append(metrics[m])
  return res

def plot_one_case(name, datasets, methods_dict, order, metric):
  try:
    output = "extreme_case_" + name + ".svg"
    ylabel = "Average unrooted RF distance"
    method_to_values = get_method_to_values(datasets, methods_dict, metric)
    plotter = boxplot.BoxPlot(name, ylabel, order)
    for m in methods_dict:
      plotter.add_elem(methods_dict[m], method_to_values[m])

    plotter.plot(output)
  except:
    print("Failed to plot " + name)

def plot_extreme_cases():
  datasets = []
  methods = []
  subst_model = "GTR+G"
  methods.append(("SpeciesRax", "speciesrax-dtl-raxml-perfam-hybrid"))
  methods.append(("MiniNJ", "njrax-mininj-raxml-ng"))
  methods.append(("AstralPro", "astralpro-raxml-ng"))
  methods.append(("CherryFusion", "njrax-cherry-raxml-ng"))
  methods.append(("NJst", "njst-original"))
  methods.append(("DupTree", "duptree"))
  methods_dict = {}
  order = []
  for t in methods:
    full_name = t[1] + "." + subst_model.lower()
    methods_dict[full_name] = t[0]
    order.append(t[0])
  metric = "species_unrooted_rf"
  for dataset_name in extreme_datasets:
    all_replicates = get_dataset_list(extreme_datasets[dataset_name], plot_replicates)
    plot_one_case(dataset_name, all_replicates, methods_dict, order, metric)

if (do_run):
  run_extreme_cases()
if (do_plot):
  plot_extreme_cases()



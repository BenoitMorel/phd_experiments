import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import run_all_species
from run_all_species import SpeciesRunFilter
import plot_speciesrax
import simulations_common
import plot_simulations

do_run = True
do_plot = not do_run
datasets = []
cores = 40
subst_model = "GTR+G"
#subst_model = "F84"
gene_trees = ["raxml-ng"]
launch_mode = "normald"
replicates = range(3000, 3050)
tag = "withoutils"
fixed_point = "ssim_" + tag + "_s50_f100_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed20"


xlabeldict = {}
xlabeldict["msmf"] = "Per species and per gene missing data"
xlabeldict["ms"] = "Per species missing data"
xlabeldict["mf"] = "Per gene missing data"
ylabel = "Average RF distance"

varying_params = []
if (do_run):
  varying_params.append((None, ["none"]))
varying_params.append(("msmf", ["ms0.2_mf0.2", "ms0.3_mf0.3", "ms0.4_mf0.4", "ms0.5_mf0.5", "ms0.6_mf0.6"]))
varying_params.append(("ms", ["ms0.4", "ms0.5", "ms0.6", "ms0.7", "ms0.8"]))
varying_params.append(("mf", ["mf0.4", "mf0.5", "mf0.6", "mf0.7", "mf0.8"])) 
#varying_params.append(("mf", ["mf0.9"]))

# NOT USED IN THE PAPER:
#varying_params.append(("families", ["f50"]))
#varying_params.append(("species", ["s25", "s100"]))
#varying_params.append(("bl", ["bl0.01", "bl0.1",  "bl10.0", "bl100.0", "bl200.0"]))
#varying_params.append(("sites", ["sites50", "sites200"]))


# metric to plot
metric_names = ["species_unrooted_rf"]
#metric_names = ["species_grf"]

# methods to plot
methods_tuples = []
methods_tuples.append(("astrid-fastme_raxml-ng", "ASTRID"))
methods_tuples.append(("astralmp_raxml-ng", "Astral"))
methods_tuples.append(("aster_raxml-ng", "Aster"))
methods_tuples.append(("fastrfs-raxml-ng_single", "FastRFS"))
methods_tuples.append(("asteroid-raxml-ng", "Asteroid"))

# run run_filter on all datasets in dataset
def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = fam.get_datadir(dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def run_varying_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()  
#  run_filter.generate = True
  run_filter.pargenes = True
  run_filter.starting_gene_trees = gene_trees
  run_filter.cleanup = True
  run_filter.asteroid = True
  run_filter.astral_mp = True
  run_filter.fastrfs = True
  run_filter.astrid_single = True
  run_filter.disable_all()  
  run_filter.asteroid = True 
  run_filter.analyze = True 
  
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates)
    run_species_methods(datasets, subst_model, cores, run_filter, launch_mode)

def plot_varying_experiment():
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates, True)
    print("Plotting parameter " + entry[0])
    for metric in metric_names:
      param = entry[0]
      output = simulations_common.get_plot_name("asteroid_" + tag, param, subst_model, metric)  
      plot_simulations.plot_varying_params(datasets, param, metric, methods_tuples, subst_model, output, xlabel = xlabeldict[param], ylabel = ylabel)


if (do_run):
  run_varying_experiment()
if (do_plot):
  plot_varying_experiment()



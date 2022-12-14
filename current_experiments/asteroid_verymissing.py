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

do_run = False
do_plot = not do_run
datasets = []
cores = 40
subst_model = "GTR+G"
gene_trees = ["raxml-ng"]
launch_mode = "normald"
replicates = range(3000, 3050)
tag = "veryhighmiss"
fixed_point = "ssim_" + tag + "_s50_f1000_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop50000000_ms0.6_mf0.6_seed20"



# pop10         0
# pop50000000  0.22
# pop100000000 0.33
# pop500000000  0.68
# pop1000000000 0.82
#
#

varying_params = []
if (do_run):
  varying_params.append((None, ["none"]))
varying_params.append(("msmf", ["ms0.5_mf0.5", "ms0.55_mf0.55",  "ms0.65_mf0.65", "ms0.7_mf0.7", "ms0.75_mf0.75"]))
varying_params.append(("population", ["pop10", "pop100000000", "pop500000000", "pop1000000000"]))
varying_params.append(("discordance", ["pop10", "pop100000000", "pop500000000", "pop1000000000"]))
varying_params.append(("families", ["f250", "f500", "f2000"]))
varying_params.append(("bl", ["bl0.05", "bl0.1",  "bl10.0", "bl100.0", "bl200.0"]))
varying_params.append(("sites", ["sites50", "sites200", "sites500"]))
varying_params.append(("species", ["s25", "s75", "s100", "s125", "s150"]))

xlabeldict = {}
xlabeldict["msmf"] = "Per species and per gene missing data"
xlabeldict["population"] = "Population size"
xlabeldict["discordance"] = "Level of discordance"
xlabeldict["families"] = "Gene number"
xlabeldict["bl"] = "Branch length scaler"
xlabeldict["sites"] = "Average sequence length"
xlabeldict["species"] = "Species number"
ylabel = "Average RF distance"


# metric to plot
metric_names = ["species_unrooted_rf"]

# methods to plot
methods_tuples = []
methods_tuples.append(("astrid-fastme_raxml-ng", "ASTRID"))
methods_tuples.append(("astralmp_raxml-ng", "ASTRAL-III"))
methods_tuples.append(("aster_raxml-ng", "ASTER"))
methods_tuples.append(("fastrfs-raxml-ng_single", "FastRFS"))
methods_tuples.append(("asteroid-raxml-ng", "Asteroid"))

# run run_filter on all datasets in dataset
def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = fam.get_datadir(dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def run_varying_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.generate = True
  run_filter.pargenes = True
  run_filter.starting_gene_trees = gene_trees
  run_filter.cleanup = True
  run_filter.asteroid = True
  run_filter.fastrfs = True
  run_filter.astrid_single = True
  run_filter.astral_mp = True
  run_filter.disable_all()  
  run_filter.analyze = True 
  run_filter.asteroid = True 
  
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates)
    run_species_methods(datasets, subst_model, cores, run_filter, launch_mode)

def plot_varying_experiment():
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list_list(fixed_point, entry[1], replicates, True)
    print("Plotting parameter " + entry[0])
    for metric in metric_names:
      param = entry[0]
      output = simulations_common.get_plot_name("asteroid_" + tag, param, subst_model, metric)  
      plot_simulations.plot_varying_params(datasets, param, metric, methods_tuples, subst_model, output, xlabel = xlabeldict[param], ylabel = ylabel)


if (do_run):
  run_varying_experiment()
if (do_plot):
  plot_varying_experiment()




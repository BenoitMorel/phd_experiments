
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
gene_trees = ["raxml-ng"]
launch_mode = "normald"
replicates = range(3000, 3100)
varying_params = []



varying_params.append((None, ["none"]))
varying_params.append(("mu", ["mu0.25", "mu0.75", "mu1.0"]))
#varying_params.append(("families", ["f200", "f1000"]))
#varying_params.append(("species", ["s15", "s35", "s50"]))

tag = "highilsmissing"
fixed_point = "ssim_" + tag + "_s25_f500_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop1000000000_mu0.5_theta5.0_seed20"

# metric to plot
metric_names = ["species_unrooted_rf"]

# methods to plot
methods_tuples = []
methods_tuples.append(("minibme-mininj_raxml-ng", "MiniBME"))
methods_tuples.append(("minibmepruned-mininj_raxml-ng", "MiniBMEPruned"))
methods_tuples.append(("njrax-mininj_raxml-ng", "MiniNJ"))
methods_tuples.append(("astral", "Astral"))
#methods_tuples.append(("fastmulrfs-single_raxml-ng", "FastMulRFS"))
#methods_tuples.append(("duptree_raxml-ng", "DupTree"))
#methods_tuples.append(("njrax-ustar_raxml-ng", "USTAR"))

# run run_filter on all datasets in dataset
def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = fam.get_datadir(dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def run_varying_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()  
  #run_filter.generate = True
  run_filter.pargenes = True
  run_filter.pargenes_starting_trees = 1
  run_filter.pargenes_bootstrap_trees = 0
  run_filter.starting_gene_trees = gene_trees
  run_filter.njrax = True
  run_filter.astrid = True
  run_filter.astral = True
  run_filter.njst = True
  run_filter.cleanup = True
  run_filter.minibme = True
  run_filter.concatenation_max = True
  run_filter.disable_all()  
  run_filter.minibmepruned = True
  run_filter.analyse = True 
  
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates)
    run_species_methods(datasets, subst_model, cores, run_filter, launch_mode)

def plot_varying_experiment():
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates, True)
    print("Plotting parameter " + entry[0])
    for metric in metric_names:
      param = entry[0]
      output = simulations_common.get_plot_name("varydl", param, subst_model, metric)  
      plot_simulations.plot_varying_params(datasets, param, metric, methods_tuples, subst_model, output)


if (do_run):
  run_varying_experiment()
if (do_plot):
  plot_varying_experiment()


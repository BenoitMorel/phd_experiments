import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import fam_data
import run_all_species
from run_all_species import SpeciesRunFilter
import plot_speciesrax

do_run = True
do_plot = False
datasets = []
cores = 40
subst_model = "GTR+G"
gene_trees = ["raxml-ng"]
launch_mode = "normald"
replicates = range(3000, 3050)
varying_params = []
varying_params += ["none"]
varying_params += ["s15", "s35", "s50"]
varying_params += ["d0.5_l0.5", "d3.0_l3.0"]
varying_params += ["t2.0", "t3.0"]
varying_params += ["sites50", "sites200", "sites500"]
varying_params += ["f50", "f300", "f1000"]

tag = "varydtl"
fixed_point = "ssim_varydtl_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop10_mu1.0_theta0.0_seed20"

# metric to plot
metric_names = ["species_unrooted_rf"]

# methods to plot
methods_tuples = []
methods_tuples.append(("speciesrax-dtl-raxml-perfam-hybrid", "SpeciesRax"))
methods_tuples.append(("astralpro-raxml-ng", "Astral-Pro"))
methods_tuples.append(("njrax-mininj-raxml-ng", "MiniNJ"))
#methods_tuples.append(("njrax-wmininj-raxml-ng", "WMiniNJ"))
methods_tuples.append(("njrax-cherry-raxml-ng", "CherryMerging"))
#methods_tuples.append(("njrax-njst-raxml-ng", "NJst"))
#methods_tuples.append(("njrax-ustar-raxml-ng", "USTAR-NJ"))
methods_tuples.append(("njrax-cherrypro-raxml-ng", "CherryMergingPro"))
#methods_tuples.append(("duptree", "DupTree"))
  
params_to_plot = ["species", "sites", "dup_rate", "families", "transfer_rate"]
fixed_params_values = {}
fixed_params_values["species"] = "30"
fixed_params_values["families"] = "100"
fixed_params_values["sites"] = "200"
fixed_params_values["tag"] = tag
fixed_params_values["dup_rate"] = "1.0"
fixed_params_values["transfer_rate"] = "1.0"


methods = []
methods_dict = {}
for t in methods_tuples:
  methods.append(t[0])
  methods_dict[t[0]] = (t[1], None)



# run run_filter on all datasets in dataset
def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = fam.get_datadir(dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

# get a list of replicated datasets with varying parameters
# from a reference dataset, a list of parameters to vary, 
# and a range of replicates
def get_dataset_list(ref_dataset, strings_to_replace, replicates):
  seeds_to_replace = []
  for i in replicates:
    seeds_to_replace.append("seed" + str(i))
  unreplicated_datasets = []
  fam_data.get_dataset_variations(unreplicated_datasets, ref_dataset, strings_to_replace)
  replicated_datasets = []
  for d in unreplicated_datasets:
    fam_data.get_dataset_variations(replicated_datasets, d, seeds_to_replace)
  return replicated_datasets

def run_varying_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.generate = True
  run_filter.starting_gene_trees = gene_trees
  run_filter.pargenes = True
  run_filter.pargenes_starting_trees = 1
  run_filter.pargenes_bootstrap_trees = 0
  run_filter.duptree = True
  run_filter.njrax = True
  run_filter.astralpro = True
  run_filter.njst = True
  run_filter.cherrypro = True
  run_filter.cherry = True
  #run_filter.speciesraxprune = True
  #run_filter.speciesraxperfamily = True
  run_filter.njrax = True
  run_filter.fastmulrfs = True 
  #run_filter.disable_all()
  run_filter.speciesraxbench = True
  #run_filter.speciesraxperfamily = True
  #run_filter.speciesraxprune = True
  run_filter.verbose = True
  # mrbayes!!
  if (False):
    run_filter.mrbayes = True
    run_filter.mb_runs = 2
    run_filter.mb_chains = 4 
    run_filter.mb_frequencies =  5000
    run_filter.mb_generations = 500000
    mb_trees = run_filter.mb_generations * run_filter.mb_runs / (run_filter.mb_frequencies)
    run_filter.mb_burnin = mb_trees / 10
  run_filter.cleanup = True
  
  datasets = get_dataset_list(fixed_point, varying_params, replicates)
  run_species_methods(datasets, subst_model, cores, run_filter, launch_mode)

def plot_varying_experiment():
  datasets = plot_speciesrax.get_datasets("ssim_" + tag)
  plot_speciesrax.generate_plot(datasets, params_to_plot, metric_names, methods, methods_dict, tag, fixed_params_values, subst_model)


if (do_run):
  run_varying_experiment()
if (do_plot):
  plot_varying_experiment()


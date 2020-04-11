import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import fam_data
import saved_metrics
import run_all_species
from run_all_species import SpeciesRunFilter

datasets = []
cores = 40
launch_mode = "normal"

do_run = False
do_plot = True

replicates = range(1000, 1020)
extreme_datasets = {}
for sites in ["100", "300"]:
# nothing 
  extreme_datasets["nothing_" + sites] = "ssim_var_s20_f100_sites" + sites + "_GTR_bl1.0_d0.0_l0.0_t0.0_p0.0_pop10_mu1.0_theta0.0_seed1000"
# DL 
  extreme_datasets["DL_" + sites] ="ssim_var_s20_f100_sites" + sites + "_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop10_mu1.0_theta0.0_seed1000"
# DTL
  extreme_datasets["DTL_" + sites] = "ssim_var_s20_f100_sites" + sites + "_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop10_mu1.0_theta0.0_seed1000"
# ILS
  extreme_datasets["ILS_" + sites] = "ssim_var_s20_f100_sites" + sites + "_GTR_bl1.0_d0.0_l0.0_t0.0_p0.0_pop100000000_mu1.0_theta0.0_seed1000"
# ILS and DL
  extreme_datasets["ILS_DL_" + sites] = "ssim_var_s20_f100_sites" + sites + "_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop100000000_mu1.0_theta0.0_seed1000"
# ILS and DTL
  extreme_datasets["ILS_DTL_" + sites] = "ssim_var_s20_f100_sites" + sites + "_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop100000000_mu1.0_theta0.0_seed1000"



def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

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

def run_extreme_cases():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.generate = True
  run_filter.pargenes_starting_trees = 10
  run_filter.pargenes_bootstrap_trees = 2
  run_filter.pargenes = True
  run_filter.duptree = True
  run_filter.njrax = True
  run_filter.astralpro = True
  run_filter.astralpromultiple = True
  run_filter.njst = True
  run_filter.concatenation_min = True
  run_filter.concatenation_max = True
  run_filter.cherry = True
  run_filter.speciesrax = True
  run_filter.speciesraxperfamily = True
  datasets = []
  for dataset_name in extreme_datasets:
    dataset = extreme_enabled[dataset]
    datasets += get_dataset_list(dataset, params, replicates)
  run_species_methods(datasets, "GTR+G", cores, run_filter, launch_mode)

def plot_extreme_cases():
  for dataset_name in extreme_datasets:
    print(dataset_name + " " + extreme_datasets[dataset_name])

if (do_run):
  run_extreme_cases()
if (do_plot):
  plot_extreme_cases()



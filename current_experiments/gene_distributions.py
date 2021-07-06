import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import run_all_species
from run_all_species import SpeciesRunFilter
from run_all import RunFilter
import plot_speciesrax
import simulations_common
import plot_simulations

datasets = []
cores = 38
subst_model = "GTR+G"
#gene_trees = ["parsimony100"]#, "parsimony500"]# "mrbayes-r2-c2-g100K-f100-b100"]
gene_trees = ["raxml-ng100", "bootstrap100", "mrbayes-r2-c2-g1M-f100-b10"]
#gene_trees = ["mrbayes-r2-c2-g1M-f1K-b100"]
#gene_trees = ["bootstrap100", "bootstrap1000"]
launch_mode = "normald"
replicates = range(3005, 3100)
varying_params = []

varying_params.append((None, ["none"]))
#varying_params.append(("sites", ["sites50", "sites200", "sites300"]))
#varying_params.append(("dup_rate", ["d0.0_l0.0", "d0.5_l0.5", "d2.0_l2.0", "d3.0_l3.0"]))
#varying_params.append(("species", ["s15", "s35", "s50", "s75"]))
#varying_params.append(("families", ["f50", "f200", "f500", "f1000"]))
#varying_params.append(("population", ["pop10000000", "pop100000000", "pop1000000000"]))

#tag = "dlsim"
#fixed_point = "ssim_" + tag + "_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop10_mu1.0_theta0.0_seed3000"
tag = "dtlsim"
fixed_point = "ssim_" + tag + "_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop10_mu1.0_theta0.0_seed3000"

generate_gene_trees = True
generate_species_trees = True

# run run_filter on all datasets in dataset
def run_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = fam.get_datadir(dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores)

def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = fam.get_datadir(dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores)


def run_varying_experiment():
  run_filter = RunFilter()
  run_filter.mrbayes = True
  run_filter.mb_frequencies = 1000
  run_filter.mb_generations = 1000000
  run_filter.mb_runs = 2
  run_filter.mb_chains = 2
  run_filter.mb_burnin = 100
  run_filter.rm_mrbayes = True
  run_filter.disable_all()
  run_filter.bootstrap = [100]
  run_filter.raxml_multi = [100]
  #run_filter.parsimony = [100, 500]
  #run_filter.bootstrap = [100]
  species_run_filter = SpeciesRunFilter()
  species_run_filter.disable_all()
  species_run_filter.starting_gene_trees = gene_trees
  species_run_filter.genetegratorbench = True
  #species_run_filter.astralpro = True
  #species_run_filter.fastmulrfs = True
  #species_run_filter.duptree = True
  #species_run_filter.cherry = True
  #species_run_filter.njrax = True
  #species_run_filter.speciesraxbench = True
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates)
    for dataset in datasets:
      dataset_dir = fam.get_datadir(dataset)
      if (generate_gene_trees):
        run_filter.run_reference_methods(dataset_dir, subst_model, cores)
      if (generate_species_trees):
        species_run_filter.run_reference_methods(dataset_dir, subst_model, cores)

def plot_varying_experiment():
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates, True)
    print("Plotting parameter " + entry[0])
    for metric in metric_names:
      param = entry[0]
      output = simulations_common.get_plot_name("varydl", param, subst_model, metric)  
      plot_simulations.plot_varying_params(datasets, param, metric, methods_tuples, subst_model, output)


run_varying_experiment()



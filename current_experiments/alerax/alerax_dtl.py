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
cores = 10
subst_model = "GTR+G"
mrbayes_trees = "mrbayes-r1-c2-g100K-f100-braxml"
gene_trees = "raxml-ng"
launch_mode = "normald"

tag = "dtl"
fixed_point = "ssim_" + tag + "_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed20"
replicates = range(3000,3050)

varying_params = []
if (True):
  #if (do_run):
  #  varying_params.append((None, ["none"]))
  #varying_params.append(("dup_rate", ["d0.5_l0.5_t0.5", "d2.0_l2.0_t2.0", "d3.0_l3.0_t3.0"]))
  #varying_params.append(("transfer_rate", ["t0.5", "t2.0", "t3.0"]))
  #varying_params.append(("population", ["pop10000000", "pop100000000", "pop1000000000"]))
  varying_params.append(("species", ["s15", "s35", "s50", "s75"]))
  #varying_params.append(("families", ["f50", "f200", "f500", "f1000"]))
  #varying_params.append(("bl", ["bl0.01", "bl0.1", "bl10.0", "bl100.0", "bl1000.0", "bl10000.0", "bl100000.0"]))
  varying_params.append(("sites", ["sites51", "sites200", "sites300"]))





methods_tuples = []
metric_names = []


if (do_plot):
  metric_names = ["species_unrooted_rf"]
  #methods_tuples.append(("astralpro2_raxml-ng","APro2-raxml"))
  methods_tuples.append(("speciesrax-dtl-raxml-ng-perfam-hybrid","SpeciesRax"))
  methods_tuples.append(("genetegrator-mininj-orunif_mrbayes-r1-c2-g100k-f100-braxml", "AleRax"))



# run run_filter on all datasets in dataset
def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = fam.get_datadir(dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def run_varying_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.starting_gene_trees = gene_trees 
  run_filter.disable_all()
  if (True):
    
    run_filter.generate = False
    run_filter.ale_undated = [mrbayes_trees]
    run_filter.alerax = [mrbayes_trees]
    #run_filter.generax_undated = True
    run_filter.analyze_gene_trees = True 
    #run_filter.fasttree = True
    #run_filter.pargenes = True
    #run_filter.plausiblerax = True
    #run_filter.mrbayes = True
    #run_filter.mb_frequencies = 100
    #run_filter.mb_generations = 100000
    #run_filter.mb_runs = 1
    #run_filter.mb_chains = 2
    #run_filter.mb_burnin = "raxml"
    #run_filter.speciesrax = [gene_trees]
    #run_filter.genetegratorbench = [mrbayes_trees]
    #run_filter.astralpro2 = [gene_trees]
  #run_filter.alegenerax_softdated = [mrbayes_trees]
  run_filter.analyze = True
  
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates)
    run_species_methods(datasets, subst_model, cores, run_filter, launch_mode)

def plot_varying_experiment():
  for entry in varying_params:
    datasets = simulations_common.get_dataset_list_list(fixed_point, entry[1], replicates, True)
    print("Plotting parameter " + entry[0])
    for metric in metric_names:
      param = entry[0]
      logscale = False
      if (metric == "runtimes" or metric == "seqtimes"):
        logscale = True
      output = simulations_common.get_plot_name("varydtl", param, subst_model, metric)  
      plot_simulations.plot_varying_params(datasets, param, metric, methods_tuples, subst_model, output,  logscale = logscale)


if (do_run):
  run_varying_experiment()
if (do_plot):
  plot_varying_experiment()





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
cores = 30
subst_model = "GTR+G"
gene_trees = ["raxml-ng"]
launch_mode = "normald"
replicates = range(3000, 3100)
varying_params = []



#varying_params.append((None, ["none"]))
#varying_params.append(("sites", ["sites50", "sites200", "sites300"]))
#varying_params.append(("families", ["f50", "f200", "f500", "f1000"]))
#varying_params.append(("bl", ["bl0.01", "bl0.1", "bl10.0", "bl100.0", "bl1000.0", "bl10000.0"]))
#varying_params.append(("dup_rate", ["d0.5_l0.5_t0.5", "d2.0_l2.0_t2.0", "d3.0_l3.0_t3.0","d5.0_l5.0_t5.0"]))
#varying_params.append(("transfer_rate", ["t0.5", "t2.0", "t3.0", "t5.0", "t10.0"]))
#varying_params.append(("dt_ratio", ["d0.0_l0.0", "d0.5_l0.5", "d2.0_l2.0", "d5.0_l5.0"]))
#varying_params.append(("population", ["pop10000000", "pop100000000", "pop1000000000"]))
#varying_params.append(("gene_conversion_rate", ["gc0.5", "gc1.0", "gc2.5", "gc10.0"]))
#varying_params.append(("species", ["s15", "s35", "s50", "s75"]))
varying_params.append(("bl", ["bl100000.0"]))

tag = "dtlsim"
fixed_point = "ssim_" + tag + "_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t1.0_gc0.0_p0.0_pop10_mu1.0_theta0.0_seed20"

# metric to plot
metric_names = ["species_unrooted_rf"]

# methods to plot
methods_tuples = []
methods_tuples.append(("minibme-mininj_raxml-ng", "MiniBME"))
methods_tuples.append(("minibmepruned-mininj_raxml-ng", "MiniBMEPruned"))
methods_tuples.append(("generax-mininj-fam_raxml-ng", "SpeciesRax"))
#methods_tuples.append(("generax-mininj-fam-fixed_raxml-ng", "SpeciesRaxFixed"))
methods_tuples.append(("njrax-mininj_raxml-ng", "MiniNJ"))
methods_tuples.append(("astralpro_raxml-ng", "Astral-Pro"))
#methods_tuples.append(("fastmulrfs-single_raxml-ng", "FastMulRFS"))
#methods_tuples.append(("duptree_raxml-ng", "DupTree"))


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

def run_varying_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.generate = True
  # mrbayes!!
  if (False):
    run_filter.mrbayes = True
    run_filter.mb_runs = 2
    run_filter.mb_chains = 4 
    run_filter.mb_frequencies =  5000
    run_filter.mb_generations = 500000
    mb_trees = run_filter.mb_generations * run_filter.mb_runs / (run_filter.mb_frequencies)
    run_filter.mb_burnin = mb_trees / 10
  run_filter.pargenes = True
  run_filter.pargenes_starting_trees = 1
  run_filter.pargenes_bootstrap_trees = 0
  run_filter.starting_gene_trees = gene_trees
  run_filter.duptree = True
  run_filter.njrax = True
  run_filter.astralpro = True
  run_filter.njst = True
  run_filter.fastmulrfs = True 
  run_filter.cleanup = True
  run_filter.speciesraxbench = True
  run_filter.minibme = True
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
      output = simulations_common.get_plot_name("varydtl", param, subst_model, metric)  
      plot_simulations.plot_varying_params(datasets, param, metric, methods_tuples, subst_model, output)


if (do_run):
  run_varying_experiment()
if (do_plot):
  plot_varying_experiment()


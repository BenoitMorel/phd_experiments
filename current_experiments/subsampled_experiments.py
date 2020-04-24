import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import fam_data
import saved_metrics
import run_all_species
import generate_families_with_subsampling
from run_all_species import SpeciesRunFilter
import plot_speciesrax

do_run = True
do_plot = False

plot_mrbayes = False
launch_mode = "normald"
cores = 40
gene_trees = ["raxml-ng"]
replicates = range(1000, 1050)
experiments = []
subsample_range = [0.01, 0.02, 0.05, 0.1]
experiments.append(("ensembl_98_ncrna_primates", "GTR+G", subsample_range))
experiments.append(("cyano_empirical", "LG+G+I", subsample_range))
plot_metric_names = ["species_unrooted_rf"]

# methods to plot
methods_tuples = []
if (not plot_mrbayes):
  methods_tuples.append(("speciesrax-dtl-raxml-perfam-hybrid", "SpeciesRax", None))
  methods_tuples.append(("astralpro-raxml-ng", "Astral-Pro", None))
  methods_tuples.append(("njrax-mininj-raxml-ng", "MiniNJ", None))
  methods_tuples.append(("njrax-cherry-raxml-ng", "CherryMerging", None))
  methods_tuples.append(("duptree-raxml-ng", "DupTree", None))
else:
  methods_tuples.append(("astralpro-raxml-ng", "Astral-Pro-ML", "blue"))
  methods_tuples.append(("njrax-mininj-raxml-ng", "MiniNJ-ML", "green"))
  methods_tuples.append(("njrax-cherry-raxml-ng", "CherryMerging-ML", "red"))
  methods_tuples.append(("astralpro-mrbayes", "Astral-Pro-MCMC", "blue"))
  methods_tuples.append(("njrax-mininj-mrbayes", "MiniNJ-MCMC", "green"))
  methods_tuples.append(("njrax-cherry-mrbayes", "CherryMerging-MCMC", "red"))



methods = []
methods_dict = {}
for t in methods_tuples:
  methods.append(t[0])
  style = None
  if ("MCMC" in t[1]):
    style = "dashed"
  methods_dict[t[0]] = (t[1], style, t[2])



def get_all_replicates(initial_dataset, sampling_ratio, replicates):
  if (sampling_ratio == 1.0):
    return [initial_dataset]
  datasets = []
  for replicate in replicates:
    dataset = initial_dataset
    dataset += "_subsample" + str(sampling_ratio)
    dataset += "_rep" + str(replicate)
    datasets.append(dataset)
  return datasets

def get_dataset_dirs(datasets):
  dataset_dirs = []
  for dataset in datasets:
    dataset_dirs.append(fam.get_datadir(dataset))
  return dataset_dirs

def run_subsampled(initial_dataset, subst_model, cores, sampling_ratio, replicates, run_filter):
  initial_datadir = fam.get_datadir(initial_dataset)
  if (sampling_ratio == 1.0):
    replicates = 1
  else:
    generate_families_with_subsampling.generate(initial_datadir, sampling_ratio, replicates)
  datasets = get_all_replicates(initial_dataset, sampling_ratio, replicates)
  dataset_dirs = get_dataset_dirs(datasets)
  for dataset_dir in dataset_dirs:
    if (not os.path.isdir(dataset_dir)):
      print("Error: directory " + dataset_dir + " does not exist")
      sys.exit(1)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def get_param_subsample(key, dataset_name):
  assert(key == "subsample")
  print(dataset_name)
  if (not "subsample" in dataset_name):
    return str(1.0)
  return dataset_name.split("_")[-2].replace(key, "")

def last_path_elem(p):
  return os.path.basename(os.path.normpath(p))

def plot_subsampled(initial_dataset, subst_model, subsamples, replicates, metric_names, methods, methods_dict, prefix):
  for metric_name in metric_names:
    output = prefix + "_" + initial_dataset + "_" + subst_model + ".svg"
    grouped_datasets = {}
    for subsample in subsamples:
      datasets = get_all_replicates(initial_dataset, subsample, replicates)
      grouped_datasets[datasets[0]] = datasets
      plot_speciesrax.plot(grouped_datasets, "subsample", methods, methods_dict, subst_model, metric_name, output, get_param_subsample, title = initial_dataset)


if (do_run):
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  #run_filter.duptree = True
  #run_filter.njst = True
  #run_filter.njrax = True
  run_filter.cherry = True
  #run_filter.astralpro = True
  #run_filter.starting_gene_trees = gene_trees
  #run_filter.speciesraxperfamily = True
  #run_filter.stag = True
  #run_filter.generaxselect = True
  run_filter.cleanup = True
  for experiment in experiments:
    dataset, subst_model, dataset_subsample_range = experiment
    for subsample in dataset_subsample_range:
      run_subsampled(dataset, subst_model, cores, subsample, replicates, run_filter)

if (do_plot):
  for experiment in experiments:
    dataset, subst_model, dataset_subsample_range = experiment
    prefix = "subsample"
    if (plot_mrbayes):
      prefix += "_mrbayes"
    plot_subsampled(dataset, subst_model, subsample_range, replicates, plot_metric_names, methods, methods_dict, prefix)


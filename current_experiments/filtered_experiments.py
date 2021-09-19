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
import generate_families_with_filter as cov_filter

do_run = False
do_plot = True
launch_mode = "normald"
cores = 40

coverages = [1.0, 0.9, 0.7, 0.5, 0.0]

run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.duptree = True
run_filter.njrax = True
run_filter.cherry = True
run_filter.njst = True
run_filter.astralpro = True
run_filter.starting_gene_trees = ["raxml-ng"]
run_filter.speciesrax = True
run_filter.speciesraxprune = True
run_filter.speciesraxperfamily = True
run_filter.stag = True
run_filter.cleanup = True
    
experiments = []
dna_model = "GTR+G"
experiments.append(("ensembl_98_ncrna_primates", dna_model, coverages))
experiments.append(("ensembl_98_ncrna_lowprimates", dna_model, coverages))
experiments.append(("ensembl_98_ncrna_mammals", dna_model, coverages))
experiments.append(("ensembl_98_ncrna_vertebrates", dna_model, coverages))
experiments.append(("ensembl_98_ncrna_sauropsids", dna_model, coverages))
#experiments.append(("cyano_empirical", "LG+G"))
#experiments.append(("cyano_empirical", "LG+G+I"))

plot_metric_names = ["species_unrooted_rf"]

methods_tuples = []  
methods_tuples.append(("speciesrax-dtl-raxml-ng-perfam-hybrid", "SpeciesRax", None))
methods_tuples.append(("speciesrax-prune", "SpeciesRaxPrune", None))
methods_tuples.append(("astralpro-raxml-ng", "Astral-Pro", None))
methods_tuples.append(("njrax-mininj-raxml-ng", "MiniNJ", None))
#methods_tuples.append(("njrax-cherry-raxml-ng", "CherryMerging", None))
#methods_tuples.append(("njrax-cherrypro-raxml-ng", "CherryMergingPro", None))


methods = []
methods_dict = {}
for t in methods_tuples:
  methods.append(t[0])
  style = None
  if ("MCMC" in t[1]):
    style = "dashed"
  methods_dict[t[0]] = (t[1], style, t[2])


def get_param_coverage(key, dataset_name):
  assert(key == "coverage")
  print(dataset_name)
  if (not "mincov" in dataset_name):
    return str(0.0)
  return dataset_name.split("mincov")[-1]

def plot_coverage(initial_dataset, subst_model, coverages, metric_names, methods, methods_dict, prefix):
  for metric_name in metric_names:
    output = prefix + "_" + initial_dataset + "_" + subst_model + ".svg"
    grouped_datasets = {}
    for coverage in coverages:
      datasets = [cov_filter.get_output_dir(initial_dataset, coverage, min_sites = 0)]
      grouped_datasets[datasets[0]] = datasets
      plot_speciesrax.plot(grouped_datasets, "coverage", methods, methods_dict, subst_model, metric_name, output, get_param_coverage, title = initial_dataset)



def run_with_coverage(input_datadir, subst_model, coverage, min_sites = 0):
  
  output_dir = cov_filter.get_output_dir(input_datadir, coverage, min_sites)
  if (not os.path.isdir(output_dir)):
    cov_filter.generate(input_datadir, coverage, min_sites)
  run_filter.run_reference_methods(output_dir, subst_model, cores, launch_mode)

if (do_run):
  for dataset in datasets:
    datadir = fam.get_datadir(dataset[0])
    subst_model = dataset[1]
    coverages = dataset[2]
    for coverage in coverages:
      run_with_coverage(datadir, subst_model, coverage, min_sites = 0)


if (do_plot):
  for experiment in experiments:
    dataset, subst_model, coverages = experiment
    prefix = "coverage"
    plot_coverage(dataset, subst_model, coverages, plot_metric_names, methods, methods_dict, prefix)







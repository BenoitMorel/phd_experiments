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
launch_mode = "normald"
cores = 30


run_inputs_aa = []
run_inputs_aa.append(("raxml-ng", "LG+G"))

run_inputs_dna = []
#run_inputs_dna.append(("fasttree", "GTR"))
#run_inputs_dna.append(("raxml-ng", "GTR+G"))
run_inputs_dna.append(("mrbayes-r2-c2-g1M-f1K-b100", "GTR+G"))


run_inputs_true = [("true", "true")]

run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.mrbayes = True
run_filter.mb_frequencies = 10000
run_filter.mb_generations = 1000000
run_filter.mb_runs = 2
run_filter.mb_chains = 2
run_filter.mb_burnin = 100
run_filter.rm_mrbayes = True

#run_filter.pargenes = True
#run_filter.concatenation_min = True
#run_filter.njrax = True
run_filter.cherry = True
run_filter.astral = True
#run_filter.fastmulrfs = True
run_filter.minibme = True
run_filter.minibmepruned = True
run_filter.astrid = True 
datasets = []

datasets.append(("stam_DNA_59", run_inputs_dna))
#datasets.append(("stam_AA_94", run_inputs_aa))

# methods to plot
methods_tuples = []
methods_tuples.append(("speciesrax-dtl-raxml-perfam-hybrid", "SpeciesRax"))
methods_tuples.append(("astralpro-raxml-ng", "Astral-Pro"))
methods_tuples.append(("njrax-mininj-raxml-ng", "MiniNJ"))
methods_tuples.append(("njrax-cherry-raxml-ng", "CherryMerging"))
methods_tuples.append(("njrax-njst-raxml-ng", "NJst"))
methods_tuples.append(("duptree-raxml-ng", "DupTree"))
    

if (do_run):
  for dataset in datasets:
    datadir = fam.get_datadir(dataset[0])
    run_inputs = dataset[1]
    for run_input in run_inputs:
      run_filter.starting_gene_trees = [run_input[0]]
      subst_model = run_input[1]
      run_filter.run_reference_methods(datadir, subst_model, cores, launch_mode)
    
methods = []
methods_dict = {}
for t in methods_tuples:
  methods.append(t[0])
  methods_dict[t[0]] = t[1]

if (do_plot):
  for dataset_tuple in datasets:
    dataset = dataset_tuple[0]
    subst_model = dataset_tuple[1]
    plot_speciesrax.plot_runtimes(dataset, subst_model, methods, methods_dict)






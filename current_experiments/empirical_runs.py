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
import numpy as np

def floatstr(my_float):
  return np.format_float_positional(my_float, trim='-')

do_run = True
do_plot = False
launch_mode = "normald"
cores = 39
single = True
minbl = -0.0000011

run_inputs_aa = []
run_inputs_aa.append(("raxml-ng", "LG+G"))
#run_inputs_aa.append(("raxml-ng", "bestAA"))

run_inputs_dna = []
run_inputs_dna.append(("raxml-ng", "GTR+G"))


run_inputs_true = [("true", "true")]

run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.pargenes = True
#run_filter.concatenation_min = True

run_filter.minbl = minbl
run_filter.asteroid = True
#run_filter.njrax = True
if (single):
  run_filter.fastmulrfs = True
  run_filter.astral_mp = True
  pass
else:
  run_filter.speciesraxbench = True
  #run_filter.fastmulrfs = True
  #run_filter.astralpro = True
  #run_filter.speciesraxprune = True
run_filter.analyze = True

datasets = []

#hey = ['22166_0', '10069_0', '10069_1', '14446_0', '10124_0', '10124_1', '10245_0', '12595_0', '12267', '18985_0', '16150_0', '10301_0', '10301_1', '16203_30', '26956_0', '18165_0', '2218_0', '2218_1', '13318_1', '27888_0', '2177_0', '23063_0']

#for h in hey:
#  datasets.append(("treebase_" + h, run_inputs_dna))
#datasets.append(("treebase_12267_0" , run_inputs_dna))
datasets.append(("ensembl105_single_maxgapratio0.8", run_inputs_dna))

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
    try:
      datadir = fam.get_datadir(dataset[0])
      run_inputs = dataset[1]
      for run_input in run_inputs:
        run_filter.starting_gene_trees = [run_input[0]]
        if (minbl > 0.0):
          run_filter.starting_gene_trees = []
          run_filter.starting_gene_trees.append(run_input[0] + "-minbl" + floatstr(minbl))
        subst_model = run_input[1]
        run_filter.run_reference_methods(datadir, subst_model, cores, launch_mode)
    except:
      print("failed to treat " + dataset[0])
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






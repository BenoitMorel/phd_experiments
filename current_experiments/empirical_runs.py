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
cores = 40
single = True

run_inputs_aa = []
run_inputs_aa.append(("raxml-ng", "LG+G"))
#run_inputs_aa.append(("raxml-ng", "bestAA"))

run_inputs_dna = []
#run_inputs_dna.append(("fasttree", "GTR"))
#run_inputs_dna.append(("true", "true"))
run_inputs_dna.append(("raxml-ng", "GTR+G"))
#run_inputs_dna.append(("raxml-ng-minbl0.0000011", "GTR+G"))


run_inputs_true = [("true", "true")]

run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.pargenes = True
#run_filter.concatenation_min = True

run_filter.minibme = True
run_filter.minibmepruned = True
run_filter.njrax = True
if (single):
  run_filter.fastmulrfs = True
  run_filter.astral = True
  run_filter.astrid_single = True 
else:
  run_filter.speciesraxbench = True
  #run_filter.fastmulrfs = True
  #run_filter.astralpro = True
  #run_filter.speciesraxprune = True
run_filter.analyze = True

datasets = []


#datasets.append(("pdb_fungi60_single", run_inputs_aa))
#datasets.append(("alternative_life92", run_inputs_aa))
datasets.append(("aa_plants_single", run_inputs_aa))
datasets.append(("aa_plants_single_maxgapratio0.5", run_inputs_aa))
datasets.append(("aa_plants_single_maxgapratio0.66", run_inputs_aa))

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






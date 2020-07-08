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

run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.pargenes = True
run_filter.duptree = True
#run_filter.fasttree = True
run_filter.njrax = True
run_filter.cherry = True
#run_filter.njst = True
run_filter.astralpro = True
run_filter.starting_gene_trees = ["raxml-ng"]
run_filter.speciesrax = True
run_filter.speciesraxprune = True
run_filter.speciesraxperfamily = True
run_filter.stag = True
    
datasets = []
dna_model = "GTR+G"
aa_model = "WAG"
datasets.append(("aa_ensembl_98_ncrna_primates", aa_model))
datasets.append(("aa_ensembl_98_ncrna_lowprimates", aa_model))
datasets.append(("aa_ensembl_98_ncrna_mammals", dna_model))
#datasets.append(("ensembl_98_ncrna_primates", dna_model))
#datasets.append(("ensembl_98_ncrna_lowprimates", dna_model))
#datasets.append(("ensembl_98_ncrna_mammals", dna_model))
#datasets.append(("ensembl_98_ncrna_vertebrates", dna_model))
#datasets.append(("ensembl_98_ncrna_sauropsids", dna_model))
#datasets.append(("cyano_empirical", "LG+G"))
#datasets.append(("cyano_empirical", "LG+G+I"))

run_filter.cleanup = True

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
    subst_model = dataset[1]
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






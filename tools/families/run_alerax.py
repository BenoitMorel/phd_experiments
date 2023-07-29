


import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/generax')
import experiments as exp
import launch_alegenerax
import fam

def run(datadir, gene_trees, subst_model, cores, additional_arguments):
  species_tree = "true"
  dataset = os.path.basename(datadir)
  resultsdir = os.path.join("AleGeneRax", dataset, species_tree + "_start_" + gene_trees, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  launch_alegenerax.run(datadir, species_tree, gene_trees, subst_model, cores, additional_arguments, resultsdir)



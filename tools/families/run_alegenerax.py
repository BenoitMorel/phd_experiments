import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/generax')
import experiments as exp
import launch_alegenerax
import fam

def run(datadir, gene_trees, subst_model, transfer_constraint, cores, additional_arguments):
  strategy = "SPR"
  species_tree = "true"
  base = "alegenerax_" + transfer_constraint.lower() + "_" + gene_trees + "_run"
  resultsdir = fam.get_run_dir(datadir, subst_model, base)
  additional_arguments.append("--transfer-constraint")
  additional_arguments.append(transfer_constraint)
  exp.reset_dir(resultsdir)
  launch_alegenerax.run(datadir, subst_model, species_tree, gene_trees, cores, additional_arguments, resultsdir)



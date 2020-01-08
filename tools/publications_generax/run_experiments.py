import subprocess
import os
import sys
import re
import common
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter
import run_all_species
from run_all_species import SpeciesRunFilter

datasets = []
cores = 40


if (True):
  datasets = []
  subst_model = "GTR"
  
  datasets.append("ssim_s20_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10000_seed13")
  #datasets.append("ssim_s40_f100_sites100_dna_d0.2_l0.2_t0.0_p0.0")
  #common.generate_all_datasets(datasets)
  run_filter = RunFilter(True, False)
  run_filter.eval_joint_ll = False
  run_filter.analyze = True
  run_filter.pargenes = False
  run_filter.mb_frequencies = 1000
  run_filter.mb_generations = 100000
  #run_filter.pargenes_starting_trees = 1
  #run_filter.pargenes_bootstrap_trees = 5
  common.run_all_reference_methods(datasets, subst_model, cores, run_filter)
 



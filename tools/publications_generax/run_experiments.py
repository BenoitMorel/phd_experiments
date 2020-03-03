import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam_data
from run_all import RunFilter
import run_all_species
from run_all_species import SpeciesRunFilter

datasets = []
cores = 40


if (True):
  datasets = []
  subst_model = "GTR"
  
  datasets.append("ssim_s20_f100_sites100_GTR_bl1.0_d0.2_l0.2_t0.2_p0.0_pop10_mu1.0_theta0.0_seed10")
  #datasets.append("ssim_s40_f100_sites100_dna_d0.2_l0.2_t0.0_p0.0")
  #fam_data.generate_all_datasets(datasets)
  run_filter = RunFilter(True, False)
  run_filter.eval_joint_ll = False
  run_filter.analyze = True
  run_filter.pargenes = True
  run_filter.mb_frequencies = 1000
  run_filter.mb_generations = 100000
  #run_filter.pargenes_starting_trees = 1
  #run_filter.pargenes_bootstrap_trees = 5
  run_filter.run_all_reference_methods(datasets, subst_model, cores)
 



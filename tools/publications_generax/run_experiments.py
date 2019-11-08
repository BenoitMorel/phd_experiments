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
  subst_model = "GTR+G"
  
  #datasets.append("ssim_s40_f100_sites100_dna_d0.2_l0.2_t0.1_p0.0")
  #datasets.append("ssim_s40_f100_sites100_dna_d0.2_l0.2_t0.0_p0.0")
  datasets.append("ssim_s15_f100_sites80_dna_d0.1_l0.1_t0.0_p0.0")
  #common.generate_all_datasets(datasets)
  run_filter = RunFilter(False, True)
  run_filter.mrbayes = False
  run_filter.rm_mrbayes = False
  run_filter.ALE = False
  run_filter.pargenes = True
  run_filter.eval_joint_ll = False
  run_filter.analyze = True
  run_filter.pargenes_starting_trees = 1
  run_filter.pargenes_bootstrap_trees = 5
  common.run_all_reference_methods(datasets, subst_model, cores, run_filter)
 



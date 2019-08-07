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
  datasets.append("jsimdtl_s19_f100_sites100_dna4_bl0.05_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsimdtl_s19_f50_sites200_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s19_f50_sites200_dna4_bl0.5_d0.1_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s41_f50_sites200_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s19_f100_sites100_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  
  
  #common.generate_all_datasets(datasets)
  run_filter = RunFilter()
  run_filter.ALE = False
  run_filter.generaxrec = 5
  run_filter.analyze = True
  run_filter.eval_joint_ll = False
  common.run_all_reference_methods(datasets, subst_model, cores, run_filter)
 



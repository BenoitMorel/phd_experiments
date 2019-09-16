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
  #datasets.append("jsimdtl_s5_f100_sites100_dna4_bl0.05_d0.1_l0.2_t0.1_p0.0")
  
  #datasets.append("jsimdtl_s19_f50_sites100_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsimdtl_s27_f50_sites100_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  datasets.append("jsim_s19_f50_sites75_dna4_bl0.5_d0.3_l0.3_t0.0_p0.0")
  #datasets.append("jsim_s41_f50_sites75_dna4_bl0.5_d0.3_l0.3_t0.0_p0.0")
  
  
  #common.generate_all_datasets(datasets)
  run_filter = RunFilter()
  run_filter.mrbayes = True
  run_filter.rm_mrbayes = False
  run_filter.ALE = True
  run_filter.analyze = True
  common.run_all_reference_methods(datasets, subst_model, cores, run_filter)
 



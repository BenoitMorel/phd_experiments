import subprocess
import os
import sys
import re
import common
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter




if (True):
  subst_models = ["LG", "LG+G"]
  datasets = ["cyano_simulated"]
  cores = 512
  do_generate = 0
  run_filter = RunFilter()
  run_filter.EXA_runs = 4
  run_filter.EXA_chains = 2
  run_filter.EXA_generations = 500000
  run_filter.EXA_frequencies = 500
  run_filter.EXA_burnin = 100
  common.submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores)




if (False):
  subst_models = ["GTR+G"]
  datasets = []
  cores = 64
  do_generate = 1
  fixed_point_dl = "jsim_s19_f100_sites250_dna4_bl0.5_d0.25_l0.25_t0.0_p0.0"
  fixed_point_dtl = "jsimdtl_s19_f100_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0"
  datasets.append(fixed_point_dl)
  datasets.append(fixed_point_dtl)
  
  common.add_dataset(datasets, fixed_point_dl, ["p0.1", "p0.2", "p0.3", "p0.5", "p0.75"])
  common.add_dataset(datasets, fixed_point_dl, ["d0.01_l0.01", "d0.05_l0.05", "d0.1_l0.1", "d0.4_l0.4"])
  common.add_dataset(datasets, fixed_point_dl, ["d0.1", "d0.2", "d0.3", "d0.4"])
  common.add_dataset(datasets, fixed_point_dl, ["sites100", "sites500", "sites750"])
  common.add_dataset(datasets, fixed_point_dl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  common.add_dataset(datasets, fixed_point_dl, ["s5", "s10", "s27", "s41"])
  
  common.add_dataset(datasets, fixed_point_dtl, ["p0.1", "p0.2", "p0.3", "p0.5"])
  common.add_dataset(datasets, fixed_point_dtl, ["s5", "s10", "s12", "s16", "s27", "s41"])
  common.add_dataset(datasets, fixed_point_dtl, ["sites100", "sites500", "sites750"])
  common.add_dataset(datasets, fixed_point_dtl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  common.add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.02_t0.01", "d0.05_l0.1_t0.05", "d0.15_l0.3_t0.15", "d0.2_l0.4_t0.2"])
  common.add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.2_t0.19", "d0.05_l0.2_t0.15", "d0.15_l0.2_t0.05", "d0.19_l0.2_t0.01"])
  
  common.submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores)
  



import subprocess
import os
import sys
import re
import common

datasets = []
cores = 32
do_generate = 1

if (False):
  #datasets.append("jsim_s19_f100_sites500_dna4_bl1.0_d0.25_l0.25_t0.0")
  #datasets.append("jsim_s12_f2_sites10_dna4_bl1.0_d0.05_l0.05_t0.0")
  #common.generate_all_datasets(datasets)
  #common.run_all_reference_methods(datasets, cores)
  #common.run_all_ALE(datasets, 1)
  #common.run_all_generax(datasets)
  #common.run_all_analyzes(datasets)
  exit(0)
  

if (True):
  common.submit_multiple_experiments_haswell(["cyano_simulated"], 0, 512)

if (False):
  fixed_point_dl = "jsim_s19_f100_sites500_dna4_bl0.5_d0.25_l0.25_t0.0"
  datasets.append(fixed_point_dl)
  #common.add_dataset(datasets, fixed_point_dl, ["d0.01_l0.01", "d0.05_l0.05", "d0.1_l0.1", "d0.4_l0.4"])
  #common.add_dataset(datasets, fixed_point_dl, ["d0.1", "d0.2", "d0.3", "d0.4"])
  #common.add_dataset(datasets, fixed_point_dl, ["sites100", "sites250", "sites750", "sites1000"])
  #common.add_dataset(datasets, fixed_point_dl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  #common.add_dataset(datasets, fixed_point_dl, ["s5", "s10", "s27", "s41"])


  common.add_dataset(datasets, fixed_point_dl, ["sites1000", "sites250", "d0.3_l0.25", "d0.1_l0.25"])

if (False):
  fixed_point_dtl = "jsimdtl_s16_f100_sites500_dna4_bl0.5_d0.1_l0.2_t0.1"
  datasets.append(fixed_point_dtl)
  common.add_dataset(datasets, fixed_point_dtl, ["s5", "s10", "s12", "s19", "s27", "s41"])
  common.add_dataset(datasets, fixed_point_dtl, ["sites100", "sites250", "sites750", "sites1000"])
  common.add_dataset(datasets, fixed_point_dtl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  common.add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.02_t0.01", "d0.05_l0.1_t0.05", "d0.15_l0.3_t0.15", "d0.2_l0.4_t0.2"])
  common.add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.2_t0.19", "d0.05_l0.2_t0.15", "d0.15_l0.2_t0.05", "d0.19_l0.2_t0.01"])


common.submit_multiple_experiments_haswell(datasets, do_generate, cores)




import subprocess
import os
import sys
import re
import common
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter

datasets = []
cores = 40

if (True):
  datasets = []
  #datasets.append("jsimdtl_s12_f200_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  datasets.append("jsim_s12_f200_sites250_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  common.generate_all_datasets(datasets)
  run_filter = RunFilter()
  run_filter.EXA_chains = 2
  run_filter.EXA_runs = 4
  run_filter.EXA_frequencies = 1000
  run_filter.EXA_generations = 1000000
  run_filter.EXA_burnin = 100
  common.run_all_reference_methods(datasets, "GTR+G", cores = 40, run_filter = run_filter)



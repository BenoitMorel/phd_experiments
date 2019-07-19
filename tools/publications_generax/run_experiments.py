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
  datasets.append("jsim_s19_f100_sites250_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  datasets.append("jsim_s27_f100_sites250_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #common.generate_all_datasets(datasets)
  species_run_filter = SpeciesRunFilter()
  #species_run_filter.pargenes = False
  common.run_species_methods(datasets, subst_model, cores = cores, run_filter = species_run_filter)


if (False):
  datasets = []
  subst_model = "GTR+G+I"
  datasets.append("jsimdtl_s5_f10_sites100_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  run_filter = RunFilter()
  common.run_all_reference_methods(datasets, subst_model, cores, run_filter)
 



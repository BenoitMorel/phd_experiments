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


def run_species_methods(datasets, subst_model, cores, run_filter):
  for dataset in datasets:
    print("Run reference species methods for " + dataset)
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    save_sdtout = sys.stdout
    redirected_file = os.path.join(dataset_dir, "species_logs_run_all." + subst_model + ".txt")
    print("Redirected logs to " + redirected_file)
    sys.stdout.flush()
    sys.stdout = open(redirected_file, "w")
    run_all_species.run_reference_methods(dataset_dir, subst_model, cores, run_filter)
    sys.stdout = save_sdtout
    print("End of run_all")
    sys.stdout.flush()

if (True):
  datasets = []
  subst_model = "GTR+G"
  #datasets.append("jsim_s12_f50_sites250_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s12_f50_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s19_f50_sites250_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s19_f50_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s27_f50_sites250_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s27_f50_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s41_f50_sites250_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s41_f50_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  
  #datasets.append("jsim_s19_f150_sites75_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s19_f150_sites75_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s12_f150_sites75_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s12_f150_sites75_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s27_f150_sites75_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s27_f150_sites75_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")
  #datasets.append("jsim_s41_f150_sites75_dna4_bl0.5_d0.2_l0.2_t0.0_p0.0")
  #datasets.append("jsimdtl_s41_f150_sites75_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0")

  #datasets.append("cyano_simulated")


  #common.generate_all_datasets(datasets)
  species_run_filter = SpeciesRunFilter()
  common.run_species_methods(datasets, subst_model, cores = cores, run_filter = species_run_filter)




import subprocess
import os
import sys
import re
import common
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter
import run_all_species
from run_all_species import SpeciesRunFilter
import discordance_rate

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

datasets = []
subst_model = "GTR"
species = range(30, 31, 10)
if (True):
  seeds = range(10, 11)
  for s in species:
      for seed in seeds:
        # only dtl
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.2_l0.2_t0.2_p0.0_pop10_seed" + str(seed))
        # only ILS
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_p0.0_pop100000_seed" + str(seed))
        # only dl
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.2_l0.2_t0.0_p0.0_pop10_seed" + str(seed))
        # idtl
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop100000_seed" + str(seed))
        # idtl publication
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10_seed" + str(seed))
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10000_seed" + str(seed))
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop20000_seed" + str(seed))
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop30000_seed" + str(seed))
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop50000_seed" + str(seed))
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop100000_seed" + str(seed))
        
        datasets.append("ssim_s" + str(15) + "_f100_sites100_GTR_bl10.0_d0.0_l0.0_t0.0_p0.0_pop10000_seed" + str(seed))

common.generate_all_datasets(datasets)

if (False): # Species tree inference
  species_run_filter = SpeciesRunFilter()
  species_run_filter.disable_all()
  species_run_filter.enable_fast_methods()
  run_species_methods(datasets, subst_model, cores = cores, run_filter = species_run_filter)
    
if (True): # Gene tree inference
  run_filter = RunFilter(True, False)
  #run_filter.disable_all()
  #run_filter.mrbayes = True
  run_filter.rm_mrbayes = False
  #run_filter.ALE = True
  #run_filter.raxml = True
  #run_filter.pargenes = True
  #run_filter.pargenes_starting_trees = 5
  #run_filter.pargenes_bootstrap_trees = 5
  #run_filter.generax = True
  run_filter.analyze = True
  
  #run_filter.eval_joint_ll = False
  #run_filter.mb_frequencies = 1000
  #run_filter.mb_generations = 10000
  common.run_all_reference_methods(datasets, subst_model, cores, run_filter)
 

#species_run_filter.speciesraxfastdl = True
#species_run_filter.pargenes = False
#species_run_filter.speciesraxfastdtl = True




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

datasets = []
subst_model = "GTR"
species = range(40, 101, 10)
if (True):
  seeds = range(10, 15)
  for s in species:
      for seed in seeds:
        # only dtl
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.2_l0.2_t0.2_p0.0_pop10_seed" + str(seed))
        # only ILS
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_p0.0_pop100000_seed" + str(seed))
        # only dl
        datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.2_l0.2_t0.0_p0.0_pop10_seed" + str(seed))

#common.generate_all_datasets(datasets)
species_run_filter = SpeciesRunFilter()
species_run_filter.disable_all()
#species_run_filter.enable_fast_methods()
species_run_filter.speciesraxfastdl = True
#species_run_filter.pargenes = False
#species_run_filter.speciesraxfastdtl = True

run_species_methods(datasets, subst_model, cores = cores, run_filter = species_run_filter)




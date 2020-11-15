import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import fam_data
import saved_metrics
import run_all_species
from run_all_species import SpeciesRunFilter
import run_all
datasets = []
cores = 40
launch_mode = "normal"




test_dataset = "ssim_test_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop10_mu1.0_theta0.0_seed20"


test_subst_model = "GTR+G"
test_datasets = [test_dataset]

def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def run_test_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.generate = True
  run_filter.pargenes = True
  run_species_methods(test_datasets, test_subst_model, cores, run_filter, launch_mode)

run_test_experiment()



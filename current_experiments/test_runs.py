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




astral_dataset = "ssim_astral_s25_f1000_sites100_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop470000000_mu1.0_theta0.0_seed20"


test_subst_model = "GTR+G"
test_datasets = []
for i in range(3000, 3050):
  test_datasets.append("ssim_varydtl_s15_f100_sites200_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop10_mu1.0_theta0.0_seed" + str(i))

def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def run_test_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.mrbayes = True
  run_filter.mb_runs = 2 
  run_filter.mb_chains = 4 
  run_filter.mb_frequencies = 250
  run_filter.mb_generations = 25000
  run_filter.mb_burnin = 10
  run_filter.starting_gene_trees = ["mrbayes"]
  run_filter.njrax = True
  run_filter.njst = True
  run_filter.cherry = True
  run_filter.astralpro= True
  run_species_methods(test_datasets, test_subst_model, cores, run_filter, launch_mode)

run_test_experiment()



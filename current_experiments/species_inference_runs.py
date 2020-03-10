import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import saved_metrics
from run_all import RunFilter
import run_all_species
from run_all_species import SpeciesRunFilter
import discordance_rate

datasets = []
cores = 40
launch_mode = "normal"

varying_enabled = False
test_enabled = True


varying_subst_model = "GTR"
varying_dataset = "ssim_vary_s20_f100_sites75_GTR_bl1.0_d0.0_l0.0_t0.0_p0.0_pop1_mu1.0_theta0.0_seed20"
varying_replicates = 10

test_subst_model = "GTR"
#test_dataset = "ssim_test_s10_f10_sites75_GTR_bl1.0_d0.3_l0.3_t0.3_p0.0_pop10_mu1.0_theta0.0_seed20"
test_dataset = "ssim_test_s20_f50_sites75_GTR_bl1.0_d0.2_l0.2_t0.1_p0.0_pop10_mu1.0_theta0.0_seed20"

def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)


def get_dataset_list(ref_dataset, varying_param, param_list, replicates):
  res = []
  return res



def run_varying_experiment():
  pass

def run_test_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.generate = True
  run_filter.pargenes = True
  run_filter.duptree = True
  run_filter.njrax = True
  run_filter.astralpro = True
  run_filter.speciesraxfastdtl = True
  run_filter.concatenation_naive = True
  run_species_methods([test_dataset], test_subst_model, cores, run_filter, launch_mode)



if (varying_enabled):
  run_varying_experiment()
if (test_enabled):
  run_test_experiment()

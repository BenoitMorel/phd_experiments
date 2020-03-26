import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import fam_data
import saved_metrics
from run_all import RunFilter
import run_all_species
from run_all_species import SpeciesRunFilter
import discordance_rate

datasets = []
cores = 40
launch_mode = "normal"

varying_enabled = True
test_enabled = False



astral_dataset = "ssim_astral_s25_f1000_sites100_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop470000000_mu1.0_theta0.0_seed20"

varying_subst_model = "GTR+G"
varying_dataset = "ssim_var_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop470000000_mu1.0_theta0.0_seed20"
#varying_dataset = "ssim_var_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t1.0_p0.0_pop10_mu1.0_theta0.0_seed20"
#varying_replicates = range(20, 30)
varying_replicates = range(30, 45)
varying_params = []
varying_params += ["none"]
varying_params += ["d0.1_t0.1_l0.1", "d0.5_t0.5_l0.5", "d2.0_l2.0_t2.0", "d3.0_l3.0_t3.0"]
varying_params += ["pop10", "pop50000000","pop100000000","pop1000000000"]
varying_params += ["t0.0", "t0.5", "t2.0", "t4.0"]
varying_params += ["sites200", "sites500"]
varying_params += ["f20", "f50", "f200", "f500", "f1000"]
varying_params += ["s15", "s35"]

test_subst_model = "GTR+G"
test_datasets = []

def run_species_methods(datasets, subst_model, cores, run_filter, launch_mode):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)

def get_dataset_list(ref_dataset, strings_to_replace, replicates):
  seeds_to_replace = []
  for i in replicates:
    seeds_to_replace.append("seed" + str(i))
  unreplicated_datasets = []
  fam_data.get_dataset_variations(unreplicated_datasets, ref_dataset, strings_to_replace)
  replicated_datasets = []
  for d in unreplicated_datasets:
    fam_data.get_dataset_variations(replicated_datasets, d, seeds_to_replace)
  return replicated_datasets

def run_varying_experiment():
  run_filter = SpeciesRunFilter()
  run_filter.disable_all()
  run_filter.generate = True
  run_filter.pargenes = True
  run_filter.duptree = True
  run_filter.njrax = True
  run_filter.astralpro = True
  run_filter.njst = True
  run_filter.speciesrax = True
  run_filter.speciesraxperfamily = True
  run_filter.concatenation_min = True
  run_filter.concatenation_max = True
  run_filter.orthogenerax = True
  #run_filter.generaxselectfam = True
  run_filter.disable_all()
  run_filter.generaxselect = True
 
  run_filter.cleanup = True
  datasets = get_dataset_list(varying_dataset, varying_params, varying_replicates)
  run_species_methods(datasets, varying_subst_model, cores, run_filter, launch_mode)

def run_test_experiment():
  run_filter = SpeciesRunFilter()
  run_species_methods(test_datasets, test_subst_model, cores, run_filter, launch_mode)



if (varying_enabled):
  run_varying_experiment()
if (test_enabled):
  run_test_experiment()

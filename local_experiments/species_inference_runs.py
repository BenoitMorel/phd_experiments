import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam_data
import fam
import saved_metrics
from run_all import RunFilter
import run_all_species
from run_all_species import SpeciesRunFilter
import discordance_rate

datasets = []
cores = 40

def get_metrics_for_datasets(datasets_prefix, metric_name):
  datasets = fam_data.get_available_datasets(datasets_prefix)
  datasets_rf_dico = {}
  datasets_runtimes_dico = {}
  total = len(datasets)
  for dataset in datasets:
    dataset_dir = os.path.join(exp.families_datasets_root, dataset)
    res = saved_metrics.get_metrics(dataset_dir, metric_name)
    if (res != None):
      if (metric_name == "runtimes"):
        for key in res: 
          if ("ALE" in key):
            res[key] = str(float(res[key]) + float(res["ExaBayes"]))
      datasets_rf_dico[dataset] = res
  return datasets_rf_dico



def run_species_methods(datasets, subst_model, cores, run_filter):
  for dataset in datasets:
    dataset_dir = os.path.join("../BenoitDatasets/families", dataset)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores)

datasets = []
subst_model = "GTR"
species = range(10, 11, 10)
seeds = range(10, 11)
if (True):
  for s in species:
      for seed in seeds:
        # only dtl
        datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.2_l0.2_t0.2_p0.0_pop10_mu1.0_theta0.0_seed" + str(seed))
        # only ILS
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_p0.0_pop100000_mu1.0_theta0.0_seed" + str(seed))
        # only dl
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.2_l0.2_t0.0_p0.0_pop10_mu1.0_theta0.0_seed" + str(seed))
        # idtl
        #datasets.append("ssim_s" + str(s) + "_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop20000_mu1.0_theta0.0_seed" + str(seed))
        pass

seeds = range(20, 25)
dlrange = (0.05, 0.1, 0.25, 0.4)
dtl_range = (0.05, 0.1, 0.25, 0.4)
pop_range = (10000, 20000, 30000, 50000)
fam_range = (10, 25, 50, 200, 500)
fam_range_missing = [100, 250, 500, 750, 1000]
bl_range = (0.01, 0.1, 2.0, 5.0, 10.0, 100.0)
mu_range =  (0.4, 0.6, 0.8, 1.0)
theta_range = (0.0, 0.05, 0.1, 0.15, 0.2)
if (False):
  for dl in dlrange:
    for seed in seeds:
      datasets.append("ssim_s40_f100_sites100_GTR_bl1.0_d" + str(dl) + "_l" + str(dl) + "_t0.0_p0.0_pop10_mu1.0_theta0.0_seed" + str(seed))
  for pop in pop_range:
    for seed in seeds:
      datasets.append("ssim_s40_f100_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_p0.0_pop" + str(pop) + "_mu1.0_theta0.0_seed" + str(seed))

  for bl in bl_range:
    for seed in seeds:
      datasets.append("ssim_s40_f100_sites100_GTR_bl" + str(bl) + "_d0.1_l0.1_t0.1_p0.0_pop10000_mu1.0_theta0.0_seed" + str(seed))

  for families in fam_range:
    for seed in seeds:
      datasets.append("ssim_s40_f" + str(families) + "_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10000_mu1.0_theta0.0_seed" + str(seed))
  for dtl in dtl_range: 
    for seed in seeds:
        datasets.append("ssim_s40_f100_sites100_GTR_bl1.0_d" + str(dtl) + "_l" + str(dtl) + "_t" + str(dtl) +"_p0.0_pop10_mu1.0_theta0.0_seed" + str(seed))
if (False):
  for mu in mu_range:
    for seed in seeds:
        datasets.append("ssim_s20_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10000_mu" + str(mu) + "_theta0.0_seed" + str(seed))
  for theta in theta_range:
    for seed in seeds:
        datasets.append("ssim_s20_f100_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10000_mu0.7_theta" + str(theta) + "_seed" + str(seed))
  for families in fam_range_missing:
    for seed in seeds:
        datasets.append("ssim_s20_f" + str(families) + "_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10000_mu0.35_theta0.0_seed" + str(seed))
  for families in fam_range_missing:
    for seed in seeds:
        datasets.append("ssim_s20_f" + str(families) + "_sites100_GTR_bl1.0_d0.1_l0.1_t0.1_p0.0_pop10000_mu0.75_theta0.15_seed" + str(seed))
    
if (False):
  datasets.append("jsimdtl_s27_f50_sites100_dna4_bl1.0_d0.1_l0.1_t0.1_p0.0")
  datasets.append("jsimdtl_s27_f50_sites100_dna4_bl1.0_d0.1_l0.2_t0.2_p0.0")
  datasets.append("jsimdtl_s41_f50_sites100_dna4_bl1.0_d0.1_l0.1_t0.1_p0.0")

  

#fam_data.generate_all_datasets(datasets)
if (True): # Species tree inference
  species_run_filter = SpeciesRunFilter()
  species_run_filter.disable_all()
  #species_run_filter.enable_fast_methods()
  #species_run_filter.pargenes = True
  #species_run_filter.concatenation_naive = True
  species_run_filter.speciesraxfastdtl = True
  #species_run_filter.speciesraxprune = True
  #species_run_filter.speciesraxperfamily = True
  species_run_filter.analyze = True

  #species_run_filter.enable_fast_methods()
  run_species_methods(datasets, subst_model, cores = cores, run_filter = species_run_filter)
    




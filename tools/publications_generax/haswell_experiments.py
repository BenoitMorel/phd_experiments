import subprocess
import os
import sys
import re
import common
sys.path.insert(0, os.path.join("tools", "families"))
from run_all import RunFilter

sys.path.insert(0, 'scripts/generax')
import scaling_generax

run_cyano_simulated = False
run_simulations = False
run_empirical = False
run_scaling = True

# CYANO SIMULATED PLOTS
if (run_cyano_simulated):
  subst_models = ["LG+G+I", "WAG"] #["LG+G+I", "DAYHOFF"]
  datasets = ["cyano_simulated"]
  cores = 512
  do_generate = 0
  run_filter = RunFilter()
  run_filter.disable_all()
  run_filter.generax = True
  run_filter.analyze = True
  run_filter.EXA_runs = 2
  run_filter.EXA_chains = 4
  run_filter.EXA_generations = 1000000
  run_filter.EXA_frequencies = 1000
  run_filter.EXA_burnin = 100
  #run_filter.eval_joint_ll = True
  common.submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores, run_filter)



# PARAMETERS SIMULATED PLOTS
if (run_simulations):
  subst_models = ["GTR+G"]
  datasets = []
  cores = 64
  do_generate = 0
  run_filter = RunFilter()
  run_filter.disable_all()
  run_filter.analyze = True
  run_filter.generax = True
  run_filter.EXA_runs = 2
  run_filter.EXA_chains = 4
  run_filter.EXA_generations = 1000000
  run_filter.EXA_frequencies = 1000
  run_filter.EXA_burnin = 100
  fixed_point_dl = "jsim_s19_f100_sites250_dna4_bl0.5_d0.25_l0.25_t0.0_p0.0"
  fixed_point_dtl = "jsimdtl_s19_f100_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0"
  datasets.append(fixed_point_dl)
  datasets.append(fixed_point_dtl)
  
  common.add_dataset(datasets, fixed_point_dl, ["p0.1", "p0.2", "p0.3", "p0.5", "p0.75"])
  common.add_dataset(datasets, fixed_point_dl, ["d0.01_l0.01", "d0.05_l0.05", "d0.1_l0.1", "d0.4_l0.4"])
  common.add_dataset(datasets, fixed_point_dl, ["d0.1", "d0.2", "d0.3", "d0.4"])
  common.add_dataset(datasets, fixed_point_dl, ["sites100", "sites500", "sites750"])
  common.add_dataset(datasets, fixed_point_dl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  common.add_dataset(datasets, fixed_point_dl, ["s5", "s10", "s27", "s41"])

  common.add_dataset(datasets, fixed_point_dtl, ["p0.1", "p0.2", "p0.3", "p0.5"])
  common.add_dataset(datasets, fixed_point_dtl, ["s5", "s10", "s12", "s16", "s27", "s41"])
  common.add_dataset(datasets, fixed_point_dtl, ["sites100", "sites500", "sites750"])
  common.add_dataset(datasets, fixed_point_dtl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
  common.add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.02_t0.01", "d0.05_l0.1_t0.05", "d0.15_l0.3_t0.15", "d0.2_l0.4_t0.2"])
  common.add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.2_t0.19", "d0.05_l0.2_t0.15", "d0.15_l0.2_t0.05", "d0.19_l0.2_t0.01"])
  
  common.submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores, run_filter)
  

# EMPIRICAL PLOTS
if (run_empirical):
  subst_models_dna = ["GTR+G"]
  datasets_dna = ["ensembl_96_ncrna_primates"]
  subst_models_prot = ["LG+G"]
  datasets_prot = ["cyano_empirical"]
  cores = 512
  do_generate = 0
  run_filter = RunFilter()
  run_filter.disable_all()
  run_filter.generax = True
  run_filter.EXA_runs = 2
  run_filter.EXA_chains = 4
  run_filter.EXA_generations = 1000000
  run_filter.EXA_frequencies = 1000
  run_filter.EXA_burnin = 100
  run_filter.eval_joint_ll = True
 
  common.submit_multiple_experiments_haswell(datasets_prot, subst_models_prot, do_generate, cores, run_filter)
  common.submit_multiple_experiments_haswell(datasets_dna, subst_models_dna, do_generate, cores, run_filter)




if (run_scaling):
  dataset = "../BenoitDatasets/families/cyano_empirical"
  starting_trees = ["raxml-ng", "random"]
  models = [1] # with or without transfers
  cores_set = [4, 8, 16, 32, 64, 128, 256, 512]
  subst_model = "LG+G"
  for cores in cores_set:
    for tree in starting_trees:
      for model in models:
        scaling_generax.launch(dataset, subst_model, tree, model, "haswell", cores)



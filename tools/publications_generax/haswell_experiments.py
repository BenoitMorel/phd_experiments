import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("scripts"))
import exp
import pickle
import fam
import fam_data
from run_all import RunFilter

sys.path.insert(0, 'scripts/generax')
import scaling_generax

run_cyano_simulated = False
run_simulations = False
run_empirical = True
run_scaling = False
run_test = False

def add_dataset(datasets, fixed_point, strings_to_replace):
  for string_to_replace in strings_to_replace:
    elems_to_replace = string_to_replace.split("_")
    split = fixed_point.split("_")
    for elem in elems_to_replace:
      param_name =  re.sub("[0-9]*[\.]*[0-9]*", "", elem)
      param_value =  re.sub("[a-zA-Z]*", "", elem)
      param_position = fam_data.get_param_position(fixed_point, param_name)
      split[param_position] = param_name + str(param_value)
    dataset = "_".join(split)
    print("Add " + dataset)
    if (dataset in datasets):
      print("duplicate: " + dataset)
      exit(1)
    datasets.append(dataset)

def submit_single_experiment_haswell(dataset, subst_model, do_generate, cores, run_filter = RunFilter()):
  command = []
  command.append("python")
  command.append(os.path.join(exp.tools_root, "publications_generax", "single_experiment.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append(str(do_generate))
  command.append(str(cores))
  results_dir = os.path.join("single_experiments", dataset)
  results_dir = exp.create_result_dir(results_dir, [])
  run_filter_file = os.path.join(results_dir, "run_filter.pickle")
  command.append(run_filter_file)
  result_msg = ""
  exp.write_results_info(results_dir, result_msg)
  pickle.dump(run_filter, open(run_filter_file, "wb"))
  submit_path = os.path.join(results_dir, "sub_generax.sh")
  if (run_filter.debug):
    exp.submit(submit_path, " ".join(command), cores, "haswelld") 
  else:
    exp.submit(submit_path, " ".join(command), cores, "haswell") 

def submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores, run_filter = RunFilter()):
  for dataset in datasets:
    for subst_model in subst_models:
      submit_single_experiment_haswell(dataset, subst_model, do_generate, cores, run_filter)

# CYANO SIMULATED PLOTS
if (run_cyano_simulated):
  subst_models = ["LG+G+I"]
  datasets = ["cyano_simulated"]
  cores = 96
  do_generate = 0
  run_filter = RunFilter()
  run_filter.disable_all()
  #run_filter.mrbayes = False
  #run_filter.deleterious = False
  #run_filter.dated_ALE = True
  run_filter.rm_mrbayes = False
  run_filter.generax = True
  run_filter.analyze = True
  run_filter.eval_joint_ll = False
  submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores, run_filter)

if (run_test):
  subst_models = ["GTR+G"]
  datasets = []
  cores = 512
  do_generate = 1
  run_filter = RunFilter()
  fixed_point_dtl = "jsimdtl_s27_f1000_sites100_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0"
  datasets.append(fixed_point_dtl)
  submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores, run_filter)


# PARAMETERS SIMULATED PLOTS
if (run_simulations):
  subst_models = ["GTR+G"]
  datasets = []
  cores = 32
  do_generate = 0
  run_filter = RunFilter()
  run_filter.disable_all()
  run_filter.eccetera = True
  run_filter.mrbayes = True
  run_filter.generax = True
  run_filter.ALE = True
  run_filter.eval_joint_ll = False
  if (True):
    fixed_point_dl = "jsim_s19_f100_sites250_dna4_bl0.5_d0.25_l0.25_t0.0_p0.0"
    fixed_point_dtl = "jsimdtl_s19_f100_sites250_dna4_bl0.5_d0.1_l0.2_t0.1_p0.0"
    datasets.append(fixed_point_dl)
    datasets.append(fixed_point_dtl)
    add_dataset(datasets, fixed_point_dl, ["p0.1", "p0.2", "p0.3", "p0.5", "p0.75"])
    add_dataset(datasets, fixed_point_dl, ["d0.01_l0.01", "d0.05_l0.05", "d0.1_l0.1", "d0.4_l0.4"])
    add_dataset(datasets, fixed_point_dl, ["d0.1", "d0.2", "d0.3", "d0.4"])
    add_dataset(datasets, fixed_point_dl, ["sites100", "sites500", "sites750"])
    add_dataset(datasets, fixed_point_dl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
    add_dataset(datasets, fixed_point_dl, ["s5", "s10", "s27", "s41"])

    add_dataset(datasets, fixed_point_dtl, ["p0.1", "p0.2", "p0.3", "p0.5"])
    add_dataset(datasets, fixed_point_dtl, ["s5", "s10", "s12", "s16", "s27", "s41"])
    add_dataset(datasets, fixed_point_dtl, ["sites100", "sites500", "sites750"])
    add_dataset(datasets, fixed_point_dtl, ["bl0.01", "bl0.05", "bl0.1", "bl0.2", "bl1.0", "bl2.0"])
    add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.02_t0.01", "d0.05_l0.1_t0.05", "d0.15_l0.3_t0.15", "d0.2_l0.4_t0.2"])
    add_dataset(datasets, fixed_point_dtl, ["d0.01_l0.2_t0.19", "d0.05_l0.2_t0.15", "d0.15_l0.2_t0.05", "d0.19_l0.2_t0.01"])
  
  submit_multiple_experiments_haswell(datasets, subst_models, do_generate, cores, run_filter)
  

# EMPIRICAL PLOTS
if (run_empirical):
  subst_models_dna = ["GTR+G"]
  #datasets_dna = ["ensembl_98_ncrna_primates", "ensembl_98_ncrna_mammals", "ensembl_98_ncrna_fishes", "ensembl_98_ncrna_murinae", "ensembl_98_ncrna_sauropsids", "ensembl_98_ncrna_vertebrates"]
  datasets_dna = ["ensembl_98_ncrna_mammals"]
  #subst_models_prot = ["LG+G"]
  #datasets_prot = ["cyano_empirical"]
  #datasets_prot = ["archaea"]
  cores = 112
  do_generate = 0
  run_filter = RunFilter()
  run_filter.disable_all()
  run_filter.pargenes = True
  run_filter.pargenes_starting_trees = 5
  run_filter.pargenes_bootstrap_trees = 10
  submit_multiple_experiments_haswell(datasets_dna, subst_models_dna, do_generate, cores, run_filter)
  #submit_multiple_experiments_haswell(datasets_prot, subst_models_prot, do_generate, cores, run_filter)
  




if (run_scaling):
  dataset = "../BenoitDatasets/families/cyano_empirical"
  starting_trees = ["raxml-ng", "random"]
  models = [1] # with or without transfers
  cores_set = [4, 8, 16, 32, 64, 128, 256, 512]
  subst_model = "LG+G+I"
  for cores in cores_set:
    for tree in starting_trees:
      for model in models:
        scaling_generax.launch(dataset, subst_model, tree, model, "haswell", cores)



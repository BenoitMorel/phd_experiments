import subprocess
import os
import sys
import re
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import simulations_common
import experiments as exp
import fam

datasets = []
cores = 40
subst_model = "GTR+G"
gene_trees = ["raxml-ng"]
launch_mode = "normald"
replicates = range(3000, 3010)
varying_params = []

#varying_params.append((None, ["none"]))
#varying_params.append(("sites", ["sites50", "sites200", "sites300"]))
#varying_params.append(("dup_rate", ["d0.0_l0.0", "d0.5_l0.5", "d2.0_l2.0"]))
varying_params.append(("species", ["s15", "s35", "s50"]))
#varying_params.append(("families", ["f50", "f200", "f500", "f1000"]))
#varying_params.append(("population", ["pop10000000", "pop100000000", "pop1000000000"]))

tag = "varydl"
fixed_point = "ssim_" + tag + "_s25_f100_sites100_GTR_bl1.0_d1.0_l1.0_t0.0_p0.0_pop10_mu1.0_theta0.0_seed20"


# run run_filter on all datasets in dataset
def run_generax(datasets, subst_model, cores, launch_mode, model, radius):
  command = []
  command.append(exp.python())
  command.append(os.path.join(exp.scripts_root, "generax/launch_generax.py"))
  command.append(dataset)
  command.append(subst_model)
  command.append("SPR")
  command.append("true")
  command.append("raxml-ng")
  command.append(launch_mode)
  command.append(str(cores))
  command.append("--rec-model")
  command.append(model)
  command.append("--max-spr-radius")
  command.append(str(radius))
  command.append("--analyze")
  command.append("no")
  print("-> Running " + " ".join(command))
  subprocess.check_call(command)

for entry in varying_params:
  datasets = simulations_common.get_dataset_list(fixed_point, entry[1], replicates)
  for dataset in datasets:
    run_generax(dataset, subst_model, cores, launch_mode, "UndatedDL", 0)
    run_generax(dataset, subst_model, cores, launch_mode, "UndatedDL", 5)




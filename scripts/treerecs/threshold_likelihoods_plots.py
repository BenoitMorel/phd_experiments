import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp

def run(datadir, script_path, resultsdir):
  command = []
  command.append("python")
  command.append(script_path)
  command.append(os.path.join(datadir, "geneTrees.newick"))
  command.append(os.path.join(datadir, "speciesTree.newick"))
  command.append(os.path.join(datadir, "alignment.txt"))
  command.append(os.path.join(datadir, "mapping.txt"))
  command.append("7")
  command.append(resultsdir)
  subprocess.check_call(command)

basedir = "HBG011000"
datadir = os.path.join(exp.datasets_root, "treerecs", "phyldog_example", basedir)
script_path = os.path.join(exp.treerecs_root, "scripts", "experiments", "plot_likelihood_vs_threshold.py")
resultsdir= os.path.join(exp.results_root, "treerecs", "threshold_likelihoods_plots", basedir)

os.makedirs(resultsdir, exist_ok=True)

result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 
run(datadir, script_path, resultsdir)


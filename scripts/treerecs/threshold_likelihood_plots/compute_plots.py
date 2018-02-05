import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
import experiments as exp

def run(genedir, script_path, resultsdir):
  os.makedirs(resultsdir, exist_ok=True)
  command = []
  command.append("python")
  command.append(script_path)
  command.append(os.path.join(genedir, gene + "_raxml_support.newick"))
  command.append(os.path.join(datadir, "speciesTree.newick"))
  command.append(os.path.join(genedir, "alignment.txt"))
  command.append(os.path.join(genedir, gene + ".map"))
  command.append("7")
  command.append(resultsdir)
  subprocess.check_call(command)

datadir = os.path.join(exp.datasets_root, "phyldog_test_dataset", "treerecs_format")
gene = "HBG011000"
genedir = os.path.join(datadir, gene)
script_path = os.path.join(exp.treerecs_root, "scripts", "experiments", "plot_likelihood_vs_threshold.py")
resultsdir= os.path.join(exp.results_root, "treerecs", "threshold_likelihood_plots", gene)

result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 
run(genedir, script_path, resultsdir)


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

datadir = os.path.join(exp.datasets_root, "treerecs")

if (len(sys.argv) != 2):
  datasets = os.listdir(datadir)
  print("Syntax error: dataset required. Suggestions of datasets: ")
  print('\n'.join(datasets))
  sys.exit(0)

basedir = sys.argv[1]

datadir = os.path.join(datadir, basedir)
script_path = os.path.join(exp.treerecs_root, "scripts", "experiments", "plot_likelihood_vs_threshold.py")
resultsdir = exp.create_result_dir(os.path.join("treerecs", "threshold_likelihoods_plots", basedir))
result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 
run(datadir, script_path, resultsdir)


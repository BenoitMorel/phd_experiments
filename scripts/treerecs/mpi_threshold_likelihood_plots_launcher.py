import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp


datadir = os.path.join(exp.datasets_root, "treerecs")

if len(sys.argv) != 4:
  datasets = os.listdir(datadir)
  print("Syntax error: python threshold_likelihoods_plots.py dataset cluster cores.\n Suggestions of datasets: ")
  print('\n'.join(datasets))
  print("Cluster can be either normal, haswell or magny")
  sys.exit(0)

dataset = sys.argv[1]
cluster = sys.argv[2]
cores = int(sys.argv[3])

resultsdir = exp.create_result_dir(os.path.join("treerecs", "threshold_likelihoods_plots", dataset, cluster + "_" + str(cores), "run"))
result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 


program = os.path.join(exp.scripts_root, "treerecs", "threshold_likelihoods_plots.py")
command = "python " + program + " " + dataset + " " + str(cores) + " " + resultsdir
submit_path = os.path.join(resultsdir, "threshold_submit.sh")

print("Results will be in " + resultsdir)
exp.submit(submit_path, command, cores, cluster) 




import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp
sys.path.insert(0, os.path.join(exp.multiraxml_root, "scripts"))
import utils_multiraxml

datadir = os.path.join(exp.datasets_root, "multi-raxml", "commands")

if (len(sys.argv) != 3):
  datasets = os.listdir(datadir)
  print("Syntax : python scripts/normalrun.py command_dir threads_number")
  print("Suggestions of datasets: ")
  print('\n'.join(datasets))
  sys.exit(0)

data = sys.argv[1]
threads = int(sys.argv[2])

resultsdir = exp.create_result_dir(os.path.join("multi-raxml", "normal_run", data))
result_msg = "multi-raxml git: \n" + exp.get_git_info(exp.multiraxml_root)
exp.write_results_info(resultsdir, result_msg) 

command = os.path.join(datadir, data, "command.txt")
warning = os.path.join(datadir, data, "warning.txt")
exp.display_warning_file(warning)

utils_multiraxml.duplicateAndRun(command, threads, resultsdir, False)




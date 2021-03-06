import os
import sys
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import experiments as exp
import subprocess
import fam



datas = {}
datas["128"]  = os.path.join(exp.datasets_root, "general/128//FastaFiles/")
datas["phyldog_example"]  = os.path.join(exp.datasets_root, "general/PhyldogDataExample/FastaFiles/")
datas["ensembl"]          = os.path.join(exp.bigdatasets_root, "ensembl_8880_15/fasta_files/")
datas["simuls"]          = os.path.join(exp.bigdatasets_root, "simuls/alignments/")
datas["simuls_higher_rate"]          = os.path.join(exp.bigdatasets_root, "simuls_higher_rate/alignments/")
datas["all"] = os.path.join(exp.bigdatasets_root, "all_1kite/split_alignment")
datas["all_filtered"] = os.path.join(exp.bigdatasets_root, "all_1kite/filtered_split_alignment")
datas["muscle"] = os.path.join(exp.bigdatasets_root, "eric_tannier/vectorbase_18/MUSCLE/")

datatypes = {}
datatypes["all"] = "aa"
datatypes["all_filtered"] = "aa"

def print_help():
  print("syntax: python fromfastadir_normal.py data cluster_mode ranks use_modeltest bs_trees [additional pargenes arguments]")
  print("  possible datas: ")
  for data in datas:
    print("    " + data)
    print("    or any family dataset in BenoitDatasets/families")
  print("  possible cluster_modes: ")
  print("    " + "normal")
  print("    " + "haswell")
  print("    " + "magny")



max_args_number = 6
if ((len(sys.argv) < max_args_number)):
  print("Error! Syntax should be " )
  print_help()
  sys.exit(0)

data = sys.argv[1]
cluster_mode = sys.argv[2]
ranks = sys.argv[3]
use_modeltest = (int(sys.argv[4]) != 0)
bootstraps_number = int(sys.argv[5])



is_aa = (data in datatypes) and (datatypes[data] == "aa")
isHaswell = (cluster_mode == "haswell")
isMagny = (cluster_mode == "magny")
runner = os.path.join(exp.pargenes_root, "pargenes", "pargenes.py")
if (is_aa):
  options = os.path.join(exp.datasets_root, "pargenes", "option_files",  "raxml_global_options_aa.txt")
else:
  options = os.path.join(exp.datasets_root, "pargenes", "option_files",  "raxml_global_options.txt")
if (data == "128"):
  options = os.path.join(exp.datasets_root, "pargenes", "option_files",  "raxml_global_options_128.txt")


additional_options = ""
if (data == "all_filtered"):
  additional_options += "--modeltest-global-parameters " + os.path.join(exp.bigdatasets_root, "all_1kite", "modeltest_options.txt")

fastafiles = ""
if (data in datas):
  fastafiles = datas[data]
else:
  fastafiles = fam.get_alignments_dir(fam.get_datadir(data))
datakey = data
resultsdir = os.path.join(exp.results_root, "pargenes", data)

if (bootstraps_number != 0):
  resultsdir = os.path.join(resultsdir, "bootstraps_" + str(bootstraps_number))
else:
  resultsdir = os.path.join(resultsdir, "no_bootstraps")
if (use_modeltest):
  resultsdir += "_modeltest"


for i in range(max_args_number, len(sys.argv)):
  arg = sys.argv[i]
  resultsdir += "_" + arg  
  
resultsdir = os.path.join(resultsdir, cluster_mode + "_" + ranks)
resultsdir = os.path.join(resultsdir, "run")
resultsdir = exp.create_result_dir(resultsdir)
result_msg = "pargenes git: \n" + exp.get_git_info(exp.pargenes_root) + "\n"
result_msg += "raxml git: \n" + exp.get_git_info(os.path.join(exp.pargenes_root, "raxml-ng")) + "\n"
result_msg += "modeltest git: \n" + exp.get_git_info(os.path.join(exp.pargenes_root, "modeltest")) + "\n"
exp.write_results_info(resultsdir, result_msg) 

command = []
command.append("python3")
command.append(runner)
command.append("-a")
command.append(fastafiles)
command.append("-o")
command.append(os.path.join(resultsdir, "pargenes_run"))
command.append("-r")
command.append(options)
if (bootstraps_number > 0):
  command.append("-b")
  command.append(str(bootstraps_number))
command.append("-c")
command.append(ranks)
if (use_modeltest):
  command.append("-m")
if (is_aa):
  command.append("-d")
  command.append("aa")
for additional in additional_options.split(" "):
  if (len(additional) > 0):
    command.append(additional)

for i in range(max_args_number, len(sys.argv)):
  command.append(sys.argv[i])

if (isHaswell):
  print("executing on haswell: " + " ".join(command))
  print("")
  print("results will be in " + resultsdir)
  print("")
  submit_path = os.path.join(resultsdir, "submit.sh")
  exp.submit_haswell(submit_path, " ".join(command), ranks)
elif (isMagny):
  print("executing on haswell: " + " ".join(command))
  print("")
  print("results will be in " + resultsdir)
  print("")
  submit_path = os.path.join(resultsdir, "submit.sh")
  exp.submit_magny(submit_path, " ".join(command), ranks)
else:
  print("executing " + " ".join(command))
  subprocess.check_call(command)


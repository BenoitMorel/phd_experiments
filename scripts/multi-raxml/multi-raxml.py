import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp
import subprocess


datas = {}
datas["phyldog_example"]  = os.path.join(exp.datasets_root, "general/PhyldogDataExample/FastaFiles/")
datas["ensembl"]          = os.path.join(exp.bigdatasets_root, "ensembl_8880_15/fasta_files/")
datas["sub_ensembl_1000"] = os.path.join(exp.bigdatasets_root, "ensembl_1000_15/fasta_files/")
datas["ensembl_buggy"] = os.path.join(exp.bigdatasets_root, "ensembl_MSA_42/fasta_files/")
datas["all"] = os.path.join(exp.bigdatasets_root, "all_1kite/split_alignment")
datas["all_filtered"] = os.path.join(exp.bigdatasets_root, "all_1kite/filtered_split_alignment")
datas["muscle"] = os.path.join(exp.bigdatasets_root, "eric_tannier/vectorbase_18/MUSCLE/")

datatypes = {}
datatypes["all"] = "aa"
datatypes["all_filtered"] = "aa"

def print_help():
  print("syntax: python fromfastadir_normal.py data cluster_mode ranks use_modeltest bs_trees [additional multi-raxml arguments]")
  print("  possible datas: ")
  for data in datas:
    print("    " + data)
  print("  possible cluster_modes: ")
  print("    " + "normal")
  print("    " + "haswell")
  print("    " + "magny")



if ((len(sys.argv) < 6) or (sys.argv[1] not in datas)):
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
runner = os.path.join(exp.multiraxml_root, "multi-raxml", "multi-raxml.py")
if (is_aa):
  options = os.path.join(exp.datasets_root, "multi-raxml", "option_files",  "raxml_global_options_aa.txt")
else:
  options = os.path.join(exp.datasets_root, "multi-raxml", "option_files",  "raxml_global_options.txt")


fastafiles = datas[data]
datakey = data
resultsdir = os.path.join(exp.results_root, "multi-raxml", data)

if (bootstraps_number != 0):
  resultsdir = os.path.join(resultsdir, "bootstraps_" + str(bootstraps_number))
else:
  resultsdir = os.path.join(resultsdir, "no_bootstraps")
if (use_modeltest):
  resultsdir += "_modeltest"

max_args_number = 6

for i in range(max_args_number, len(sys.argv)):
  arg = sys.argv[i]
  resultsdir += "_" + arg  
  
resultsdir = os.path.join(resultsdir, cluster_mode + "_" + ranks)
resultsdir = os.path.join(resultsdir, "run")
resultsdir = exp.create_result_dir(resultsdir)
result_msg = "multi-raxml git: \n" + exp.get_git_info(exp.multiraxml_root) + "\n"
result_msg += "raxml git: \n" + exp.get_git_info(os.path.join(exp.multiraxml_root, "raxml-ng")) + "\n"
result_msg += "modeltest git: \n" + exp.get_git_info(os.path.join(exp.multiraxml_root, "modeltest")) + "\n"
exp.write_results_info(resultsdir, result_msg) 

command = []
command.append("python3")
command.append(runner)
command.append("-a")
command.append(fastafiles)
command.append("-o")
command.append(os.path.join(resultsdir, "multiraxml_run"))
command.append("-r")
command.append(options)
command.append("-b")
command.append(str(bootstraps_number))
command.append("-c")
command.append(ranks)
if (use_modeltest):
  command.append("-m")
if (is_aa):
  command.append("-d")
  command.append("aa")

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


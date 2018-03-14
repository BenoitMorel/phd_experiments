import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp
import subprocess


datas = {}
datas["phyldog_example"]  = os.path.join(exp.datasets_root, "general/PhyldogDataExample/FastaFiles/")
datas["ensembl"]          = os.path.join(exp.bigdatasets_root, "ensembl_8880_15/fasta_files/")
datas["sub_ensembl_1000"] = os.path.join(exp.bigdatasets_root, "ensembl_1000_15/fasta_files/")


def print_help():
  print("syntax: python fromfastadir_normal.py data cluster_mode ranks ")
  print("  possible datas: ")
  for data in datas:
    print("    " + data)
  print("  possible cluster_modes: ")
  print("    " + "normal")
  print("    " + "haswell")



if ((len(sys.argv) != 4) or (sys.argv[1] not in datas)):
  print("Error! Syntax should be " )
  print_help()
  sys.exit(0)

data = sys.argv[1]
cluster_mode = sys.argv[2]
ranks = sys.argv[3]

isHaswell = (cluster_mode == "haswell")
implementation = "--split_scheduler"
runner = os.path.join(exp.multiraxml_root, "scripts", "multiraxml_fromfastadir.py")
raxmlbin = os.path.join(exp.raxml_root, "bin")
options = os.path.join(exp.multiraxml_root, "examples", "raxml_options.txt")

fastafiles = datas[data]
resultsdir = os.path.join(exp.results_root, "multi-raxml", "fromfastadir", cluster_mode + "_" + ranks, data)
resultsdir = exp.create_result_dir(resultsdir)
result_msg = "multi-raxml git: \n" + exp.get_git_info(exp.multiraxml_root) + "\n"
result_msg += "raxml git: \n" + exp.get_git_info(exp.raxml_root) + "\n"
exp.write_results_info(resultsdir, result_msg) 

command = []
command.append("python")
command.append(runner)
command.append(implementation)
command.append(raxmlbin)
command.append(fastafiles)
command.append(resultsdir)
command.append(options)
command.append(ranks)

if (isHaswell):
  print("executing on haswell: " + " ".join(command))
  print("results will be in " + resultsdir)
  submit_path = os.path.join(resultsdir, "submit.sh")
  exp.submit_haswell(submit_path, " ".join(command), ranks)
else:
  print("executing " + " ".join(command))
  subprocess.check_call(command)


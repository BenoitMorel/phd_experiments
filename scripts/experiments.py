# Common variables used in the scripts

import os
import datetime
import sys
import subprocess

def get_parent_path(path):
  return os.path.abspath(os.path.join(path, os.pardir))

root = get_parent_path(get_parent_path(os.path.realpath(__file__)))

# scripts
scripts_root = os.path.join(root, "scripts")

# programs
programs_root = os.path.join(root, "programs")

# datasets
datasets_root = os.path.join(root, "datasets")

# results
results_root = os.path.join(root, "results")

# externals
treerecs_root = os.path.join(root, "..", "Treerecs")
treerecs_exec = os.path.join(treerecs_root, "build", "bin", "Treerecs")
multiraxml_root = os.path.join(root, "..", "multi-raxml")
multiraxml_exec = os.path.join(multiraxml_root, "build", "multi-raxml")
raxml_root = os.path.join(root, "..", "raxml-ng")

# utils
def get_git_info(repo_path):
  initial_dir = os.getcwd()
  os.chdir(repo_path)
  repo = str(subprocess.check_output(["git", "config", "--get", "remote.origin.url"]))[2:-3]
  branch = str(subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]))[2:-3]
  commit = str(subprocess.check_output(["git", "rev-parse", "HEAD"]))[2:-3]
  os.chdir(initial_dir)
  result = "Git repo: " + repo + "\nGit branch: " + branch + "\nGit commit: " + commit + "\n"
  return result

def write_results_info(resultsdir, msg):
  filename = os.path.join(resultsdir, "results.info")
  with open(filename, "w") as writer:
    writer.write("Start: "+ str(datetime.datetime.now()) + "\n")
    writer.write("Command: " + ' '.join(sys.argv) + "\n")
    writer.write("\n")
    writer.write("Experiment git: \n")
    writer.write(get_git_info("."))
    writer.write("\n")
    writer.write(msg)

def display_warning_file(warning_filename):
  if (not os.path.isfile(warning_filename)):
      return
  with open(warning_filename) as f:
    print("##########################")
    print("# Warning file content:  #")
    print("##########################")
    print("")
    print(f.read())
    print("########################")
    print("# End of warning file  #")
    print("########################")

def create_result_dir(suffix):
  base = os.path.join(results_root, suffix) + "_"
  result_dir = ""
  for i in range(0, 10000):
    result_dir = base + str(i)
    if (not os.path.isdir(result_dir)):
      os.makedirs(result_dir)
      return os.path.abspath(result_dir)

def redirect_logs(result_dir):
    logs = os.path.join(result_dir, "logs.txt")
    err = os.path.join(result_dir, "err.txt")
    print("logs redirected to " + logs)
    sys.stdout = open(logs, 'w')
    sys.stderr = open(err, 'w')

def submit_haswell(submit_file_path, command, threads):
  nodes = str((int(threads) - 1) // 16 + 1)
  with open(submit_file_path, "w") as f:
    f.write("#!/bin/bash\n")
    f.write("#SBATCH -o " + submit_file_path + ".out" + "\n")
    f.write("#SBATCH -B 2:8:1\n")
    f.write("#SBATCH -N " + str(nodes) + "\n")
    f.write("#SBATCH -n " + str(threads) + "\n")
    f.write("#SBATCH --threads-per-core=1\n")
    f.write("#SBATCH --cpus-per-task=1\n")
    f.write("#SBATCH --hint=compute_bound\n")
    f.write("#SBATCH -t 24:00:00\n")
    f.write("\n")
    f.write(command)
  subprocess.check_call(["sbatch", "-s" ,submit_file_path])



  



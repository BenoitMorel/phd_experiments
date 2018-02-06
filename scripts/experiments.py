# Common variables used in the scripts

import os
import datetime
import sys
import subprocess
root = "./"

# scripts
scripts_root = os.path.join(root, "scripts")

# datasets
datasets_root = os.path.join(root, "datasets")

# results
results_root = os.path.join(root, "results")

# externals
treerecs_root = os.path.join("..", "Treerecs")
multiraxml_root = os.path.join("..", "multi-raxml", "deps", "raxml-ng")

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


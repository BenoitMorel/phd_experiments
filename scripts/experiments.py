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

#tools
tools_root = os.path.join(root, "tools")
rf_distance_tool = os.path.join(tools_root, "treerecs", "rf_distance.py")
analyse_pargenes_jointsearch_tool = os.path.join(tools_root, "jointsearch", "analyse_pargenes_results.py")

# programs
programs_root = os.path.join(root, "programs")

# datasets
datasets_root = os.path.join(root, "datasets")
# results
results_root = os.path.join(root, "results")


# externals
github_root = os.path.join(root, "..")
benoit_datasets_root = os.path.join(github_root, "BenoitDatasets")
families_datasets_root  = os.path.join(benoit_datasets_root, "families")
treerecs_root = os.path.join(github_root, "Treerecs")
treerecs_exec = os.path.join(treerecs_root, "build", "bin", "Treerecs")
joint_likelihood_evaluator_exec = os.path.join(treerecs_root, "build", "bin", "misc", "JointLikelihoodEvaluator")
treerecs_joint_search_exec = os.path.join(treerecs_root, "build", "bin", "misc", "jointTreeSearch")
joint_search_root = os.path.join(github_root, "JointSearch")
joint_search_exec = os.path.join(joint_search_root, "build", "bin", "JointSearch")
joint_search_lib = os.path.join(joint_search_root, "build_lib", "src", "libJointSearch.so")
pargenes_root = os.path.join(github_root, "pargenes")
pargenes_script = os.path.join(pargenes_root, "pargenes", "pargenes.py")
mpischeduler_root = os.path.join(github_root, "MPIScheduler")
mpischeduler_exec = os.path.join(mpischeduler_root, "build", "mpi-scheduler")
raxml_root = os.path.join(github_root, "raxml-ng")
oldraxml_root = os.path.join(github_root, "standard-RAxML")
oldraxml_exec = os.path.join(oldraxml_root, "raxmlHPC-AVX")
bigdatasets_root = os.path.join(github_root, "datasets")
phyldog_root = os.path.join("/home/morelbt/github/PHYLDOG")
phyldog_light_exec = os.path.join(phyldog_root, "build", "bin", "phyldog_light")
zombi_script = os.path.join(github_root, "ZOMBI", "Zombi.py")
jprime_jar = os.path.join(github_root, "jprime", "jprime-0.3.6.jar")
seq_gen_exec = os.path.join(github_root, "Seq-Gen-1.3.4", "source", "seq-gen")
# constants
mpi_scheduler_heuristic = "--split-scheduler"

# utils
def get_git_info(repo_path):
  initial_dir = os.getcwd()
  os.chdir(repo_path)
  repo = str(subprocess.check_output(["git", "config", "--get", "remote.origin.url"]))[2:-3]
  branch = str(subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]))[2:-3]
  commit = str(subprocess.check_output(["git", "rev-parse", "HEAD"]))[2:-3]
  diff = "\t\t" + subprocess.check_output(["git", "diff"]).decode("ASCII").replace("\n", "\n\t\t")
  os.chdir(initial_dir)
  result = "Git repo: " + repo + "\nGit branch: " + branch + "\nGit commit: " + commit + "\n"
  result += "diff:\n " + diff + "\n"
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

def submit_haswell(submit_file_path, command, threads, debug):
  threads = int(threads)
  nodes = str((int(threads) - 1) // 16 + 1)
  logfile = os.path.join(os.path.dirname(submit_file_path), "logs.out")
  with open(submit_file_path, "w") as f:
    f.write("#!/bin/bash\n")
    f.write("#SBATCH -o " + logfile + "\n")
    f.write("#SBATCH -B 2:8:1\n")
    f.write("#SBATCH -N " + str(nodes) + "\n")
    f.write("#SBATCH -n " + str(threads) + "\n")
    f.write("#SBATCH --threads-per-core=1\n")
    f.write("#SBATCH --cpus-per-task=1\n")
    f.write("#SBATCH --hint=compute_bound\n")
    if (debug):
      f.write("#SBATCH -t 2:00:00\n")
    else:
      f.write("#SBATCH -t 24:00:00\n")

    f.write("\n")
    f.write(command)
  command = []
  command.append("sbatch")
  if (debug):
    command.append("--qos=debug")
  command.append("-s")
  command.append(submit_file_path)
  subprocess.check_call(command)

def submit_magny(submit_file_path, command, threads):
  #nodes = str((int(threads) - 1) // 16 + 1)
  #logfile = os.path.join(os.path.dirname(submit_file_path), "logs.out")
  with open(submit_file_path, "w") as f:
    f.write("#!/bin/bash\n")
    f.write("#$ -cwd -V\n")                              # Shift directories and export variables
    f.write("#$ -q bridge.q\n")            # Select the queue
    f.write("#$ -pe mvapich16 " + str(threads) + "\n")                          # Set the parallel environment
    f.write("#$ -l h_rt=24:00:00\n")                     # Request the time for the job
    f.write("#$ -N cyano_bridge\n")
    f.write("\n")
    f.write(command)
  command = []
  command.append("qsub")
  #if (int(threads) <= 32):
  #  command.append("--qos=debug")
  #command.append("-s")
  command.append(submit_file_path)
  subprocess.check_call(command)

def submit_normal(submit_file_path, command, log_cout):
    commands_list = command.split("\n")
    logfile = os.path.join(os.path.dirname(submit_file_path), "logs.out")
    for subcommand in commands_list:
      if (log_cout):
        subprocess.check_call(subcommand, shell=True)
      else:
        subprocess.check_call(subcommand + " &>> " + logfile , shell=True)


def submit(submit_file_path, command, threads, cluster):
  if (cluster == "normal"):
    submit_normal(submit_file_path, command, False)
  elif (cluster == "normald"):
    submit_normal(submit_file_path, command, True)
  elif (cluster == "haswell"):
    submit_haswell(submit_file_path, command, threads, False)
  elif (cluster == "haswelld"):
    print("debug mode")
    submit_haswell(submit_file_path, command, threads, True)
  elif (cluster == magny):
    submit_magny(submit_file_path, command, threads)
  else:
    print("unknown cluster " + cluster)
    sys.exit(1)

  




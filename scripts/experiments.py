# Common variables used in the scripts

import os
import datetime
import sys
import subprocess
import shutil

def get_parent_path(path):
  return os.path.abspath(os.path.join(path, os.pardir))

root = get_parent_path(get_parent_path(os.path.realpath(__file__)))

# scripts
scripts_root = os.path.join(root, "scripts")

#tools
tools_root = os.path.join(root, "tools")
rf_distance_tool = os.path.join(tools_root, "trees", "rf_distance.py")
rf_cells_tool = os.path.join(tools_root, "families", "rf_cells.py")
treedist_R_script = os.path.join(tools_root, "families", "treedist.R")

# programs
programs_root = os.path.join(root, "programs")

# datasets
datasets_root = os.path.join(root, "datasets")
# results
results_root = os.path.join(root, "results")
# install
installer_root = os.path.join(root, "installer")


# externals
github_root = os.path.join(root, "..")
benoit_datasets_root = os.path.join(github_root, "BenoitDatasets")
families_datasets_root  = os.path.join(benoit_datasets_root, "families")
dna4_model_samples  = os.path.join(benoit_datasets_root, "DNA_models_sample")
treerecs_root = os.path.join(github_root, "Treerecs")
treerecs_exec = os.path.join(treerecs_root, "build", "bin", "Treerecs")
joint_likelihood_evaluator_exec = os.path.join(treerecs_root, "build", "bin", "misc", "JointLikelihoodEvaluator")
treerecs_joint_search_exec = os.path.join(treerecs_root, "build", "bin", "misc", "jointTreeSearch")
treefix_exec = "treefixDTL"
joint_search_root = os.path.join(github_root, "GeneRax")
generax_exec = os.path.join(joint_search_root, "build", "bin", "generax")
speciesrax_exec = os.path.join(joint_search_root, "build", "bin", "speciesrax")
generax_gprof_exec = os.path.join(joint_search_root, "gprof_build", "bin", "generax")
generax_scalasca_exec = os.path.join(joint_search_root, "scalasca_build", "bin", "generax")
joint_search_exec = os.path.join(joint_search_root, "build", "bin", "JointSearch")
joint_search_gprof_exec = os.path.join(joint_search_root, "gprof_build", "bin", "JointSearch")
joint_search_scalasca_exec = os.path.join(joint_search_root, "scalasca_build", "bin", "JointSearch")
joint_search_lib = os.path.join(joint_search_root, "build_lib", "src", "JointSearch", "libJointSearch.so")
pargenes_root = os.path.join(github_root, "pargenes")
pargenes_script = os.path.join(pargenes_root, "pargenes", "pargenes-hpc.py")
pargenes_script_debug = os.path.join(pargenes_root, "pargenes", "pargenes-hpc-debug.py")
mpischeduler_root = os.path.join(github_root, "MPIScheduler")
mpischeduler_exec = os.path.join(mpischeduler_root, "build", "mpi-scheduler")
rfdistance_exec = os.path.join(github_root, "RaxmlRFDistance", "bin", "rfdistance-mpi")
raxml_root = os.path.join(github_root, "raxml-ng")
raxml_exec = os.path.join(raxml_root, "bin", "raxml-ng-mpi")
rf_root = os.path.join(github_root, "MPIRaxmlRFDistance")
rf_exec = os.path.join(rf_root, "bin", "rfdistance")
raxml_nompi_exec = os.path.join(raxml_root, "bin", "raxml-ng")
oldraxml_root = os.path.join(github_root, "standard-RAxML")
oldraxml_exec = os.path.join(oldraxml_root, "raxmlHPC-AVX")
bigdatasets_root = os.path.join(github_root, "datasets")
phyldog_root = os.path.join(github_root, "PHYLDOG")
phyldog_light_exec = os.path.join(phyldog_root, "build", "bin", "phyldog_light")
phyldog_exec = os.path.join(phyldog_root, "build", "bin", "phyldog")
zombi_script = os.path.join(github_root, "ZOMBI", "Zombi.py")
jprime_jar = os.path.join(github_root, "jprime", "jprime-0.3.6.jar")
jprime_deleterious_jar = os.path.join(github_root, "jprime", "jprime_0.3.5c.jar")
ale_root = os.path.join(github_root, "ALE")
ale_simul_exec = os.path.join(ale_root, "build", "bin", "simulation")
seq_gen_exec = os.path.join(github_root, "Seq-Gen-1.3.4", "source", "seq-gen")
notung_jar = os.path.join(github_root, "Notung-2.9", "Notung-2.9.jar")
mafft_exec = os.path.join(github_root, "mafft", "bin", "mafft")
stag_script = os.path.join(github_root, "STAG", "stag", "stag.py")
phylobayes_exec = os.path.join(github_root, "phylobayes4.1c", "data", "pb")
exabayes_exec = os.path.join(github_root, "exabayes-1.5", "yggdrasil")
mrbayes_exec = os.path.join(github_root, "MrBayes", "src", "mb")
ale_observe_exec = os.path.join(github_root, "ALE", "build", "bin", "ALEobserve")
ale_ml_exec = os.path.join(github_root, "ALE", "build", "bin", "ALEml_undated")
ale_ml_dated_exec = os.path.join(github_root, "ALE", "build", "bin", "ALEml")
deco_exec = os.path.join(github_root, "DeCo", "DeCo")
decostar_exec = os.path.join(github_root, "DeCoSTAR", "bin/DeCoSTAR")
guenomu_exec = os.path.join(github_root, "guenomu", "src", "guenomu")
eccetera_root = os.path.join(github_root, "ecceTERA")
eccetera_exec = os.path.join(eccetera_root, "build", "bin", "ecceTERA")
simphy_exec = os.path.join(github_root, "SimPhy_1.0.2", "bin", "simphy_lnx64")
simphy_indelible_wrapper = os.path.join(github_root, "SimPhy_1.0.2", "scripts", "INDELIble_wrapper.pl")
duptree_exec = os.path.join(github_root, "duptree", "duptree")
fastrfs_exec = os.path.join(github_root, "FastRFS", "build", "FastRFS")
astral_jar = os.path.join(github_root, "ASTRAL", "Astral", "astral.jar")
prepare_fastrfs_script = os.path.join(tools_root, "rfs", "prepareTrees.py")
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

def create_result_dir(suffix, additional_args = []):
  base = os.path.join(results_root, suffix) 
  for arg in additional_args:
    base += "_" + arg
  base += "_"
  result_dir = ""
  for i in range(0, 10000):
    result_dir = base + str(i)
    if (not os.path.isdir(result_dir)):
      os.makedirs(result_dir)
      print("Results directory: " + result_dir)
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
      print("debug on")
      f.write("#SBATCH -t 2:00:00\n")
    else:
      print("debug of")
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
    print("no debug mode")
    submit_haswell(submit_file_path, command, threads, False)
  elif (cluster == "haswelld"):
    print("debug mode")
    submit_haswell(submit_file_path, command, threads, True)
  elif (cluster == magny):
    submit_magny(submit_file_path, command, threads)
  else:
    print("unknown cluster " + cluster)
    sys.exit(1)

  
def try_make_dir(dir_name):
  try:
    os.makedirs(dir_name)
  except:
    pass

def mkdir(dir_name):
  try_make_dir(dir_name)

def relative_symlink(src, dest):
  relative_path = os.path.relpath(src, os.path.dirname(dest))
  tmp = dest + ".sym"
  os.symlink(relative_path,  tmp)
  shutil.move(tmp, dest)

def reset_dir(dir_name):
  shutil.rmtree(dir_name, True)
  os.makedirs(dir_name)

def checkAndDelete(arg, arguments):
  if (arg in arguments):
    arguments.remove(arg)
    return True
  return False

def getAndDelete(arg, arguments, default_value):
  print ("looking for " + arg + " in " + str(arguments))
  if (arg in arguments):
    index = arguments.index(arg)
    res = arguments[index + 1]
    del arguments[index + 1]
    del arguments[index]
    return res
  else:
    return default_value

def run_with_scheduler(executable, command_file, parallelization, cores, scheduler_output_dir, logname = None):
  command = ""
  out = sys.stdout
  if (logname != None):
    out = open(os.path.join(scheduler_output_dir, logname), "w")

  isMPI = (parallelization == "onecore") or (parallelization == "split")
  if (isMPI):
    command += "mpirun -np " + str(cores) + " "
  command += mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += str(cores) + " "
  command += executable + " "
  command += command_file + " "
  command += scheduler_output_dir 
  print("Running " + command)
  subprocess.check_call(command.split(" "), stdout = out, stderr = out)


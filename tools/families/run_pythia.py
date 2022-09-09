import os
import sys
import subprocess
import shutil
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import fam
import saved_metrics

def get_pythia_path(pythia_command):
  d = subprocess.check_output(["which", pythia_command])
  res = d.split("\n")[0]
  return res

def extract(datadir, run_dir):
  logs = os.path.join(run_dir, "per_job_logs")
  
  for family in fam.get_families_list(datadir):
    log = os.path.join(logs, family + "_out.txt")
    s =open(log).read()
    score = s.split()[-1][:-1]
    open(fam.get_pythia_score_path(datadir, family), "w").write(score)
    print(fam.get_pythia_score(datadir, family))

def generate_scheduler_commands_file(datadir, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      ali = fam.get_alignment(datadir, family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("-m")
      command.append(ali)
      command.append("-o")
      command.append(fam.get_pythia_score_path(datadir, family))
      command.append("-prec")
      command.append("5")
      command.append("--raxmlng")
      command.append(exp.raxml_exec_no_mpi)
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def run_pythia(datadir, cores):
  output_dir = fam.get_run_dir(datadir, "", "pythia_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, cores, output_dir)
  start = time.time()
  p = get_pythia_path(exp.pythia_exec)
  exp.run_with_scheduler(p, scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, "pythia", (time.time() - start), "runtimes") 
  extract(datadir, output_dir)


if (__name__== "__main__"):
  max_args_number = 3
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir cores.")
    sys.exit(0)


  datadir = sys.argv[1]
  cores = int(sys.argv[2])

  run_pythia(datadir, cores)




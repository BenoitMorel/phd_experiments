import subprocess
import sys
import os
sys.path.insert(0, 'scripts')
import experiments as exp

def run_scheduler(command_file, exec_path, cores, output_dir, logs_name):
  command = ""
  parallelization = "onecore"
  command += "mpirun -np " + str(cores) + " "
  command += exp.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += exec_path + " "
  command += command_file + " "
  command += output_dir + " " 
  command += "0"
  print(command.split(" "))
  logs_file = os.path.join(output_dir, logs_name)
  print("Redirecting logs to " + logs_file)
  with open(logs_file, "w") as writer:
    subprocess.check_call(command.split(" "), stdout = writer, stderr = writer)
  return command 


import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp

def print_help():
  print("syntax: python generate_parsimony.py runname sequence model trees_number")

if (len(sys.argv) != 5):
  print_help()
  sys.exit(0)

runname = sys.argv[1]
sequence = sys.argv[2]
model = sys.argv[3]
trees_number = int(sys.argv[4])

resultsdir = os.path.join("oldraxml", "generate_parsimony", model, runname)
resultsdir = exp.create_result_dir(resultsdir)

result_msg = "old raxml git: \n" + exp.get_git_info(exp.oldraxml_root) + "\n"
exp.write_results_info(resultsdir, result_msg) 

for parsimonySeed in range(1, trees_number):
  command = []
  command.append(exp.oldraxml_exec)
  command.append("-y") # stop after parsimony generation
  command.append("-m")
  command.append(model)
  command.append("-s")
  command.append(sequence)
  command.append("-p")
  command.append(str(parsimonySeed))
  command.append("-w")
  command.append(resultsdir)
  command.append("-n")
  command.append("parsi" + str(parsimonySeed))
  print("Running " + str(" ".join(command)))
  subprocess.check_call(command)


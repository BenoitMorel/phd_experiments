import subprocess
import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp

def remove(file_name):
  try:
    os.remove(file_name)
  except:
    pass

def get_tca(trees_file, run_name):
  raxml_log = "RAxML_info." + run_name
  raxml_majority_file = "RAxML_MajorityRuleExtendedConsensusTree_IC." + run_name
  redirected_log = "tca." + run_name
  generated_files = []
  generated_files.append(raxml_log)
  generated_files.append(raxml_majority_file)
  generated_files.append(redirected_log)
  #clean if needed
  for generated_file in generated_files:
    remove(generated_file)
  #prepare command
  executable = exp.oldraxml_exec
  command = []
  command.append(executable)
  command.append("-m")
  command.append("GTRCAT")
  command.append("-L")
  command.append("MRE")
  command.append("-z")
  command.append(trees_file)
  command.append("-n")
  command.append(run_name)
  #execute command
  with open(redirected_log, "w") as writer:
    subprocess.check_call(command, stdout = writer)
  # parse lgos
  raxml_lines = open(raxml_log).readlines()
  result = 0.0
  for line in raxml_lines:
    if (line.startswith("Relative tree certainty including")):
      result = float(line.split(" ")[-1][:-1])
      break
  # todo remove generated files
  for generated_file in generated_files:
    remove(generated_file)
  return result

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: python script.py trees_file run_name")
    exit(1)
  trees_file = sys.argv[1]
  run_name = sys.argv[2]
  print("TCA: " + str(get_tca(trees_file, run_name)))


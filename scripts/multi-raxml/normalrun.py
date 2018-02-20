import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp
import subprocess

# returns the new command_filename
def copyAndEditCommand(command_filename, resultsdir):
    newcommand_filename = os.path.join(resultsdir, "command.txt")
    lines = open(oldcommand_filename, "r").readlines()
    with open(newcommand_filename, "w") as writer:
        for line in lines:
            writer.write(line.replace("RESULTS_DIR", resultsdir))
    return newcommand_filename

def runCommand(command_filename, resultsdir, ranks):
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append("1")
    command.append(exp.multiraxml_exec)
    command.append(exp.multiraxml_heuristic)
    command.append(command_filename)
    command.append(resultsdir)
    command.append(str(ranks))
    print("running " + ' '.join(command))
    subprocess.check_call(command)

datadir = os.path.join(exp.datasets_root, "multi-raxml", "commands")

if (len(sys.argv) != 3):
  datasets = os.listdir(datadir)
  print("Syntax: python normalrun.py dataset ranksNumber")
  print("  Suggestions of datasets: ")
  print("  - " + '\n  - '.join(datasets))
  sys.exit(0)

basedir = sys.argv[1]
ranks = sys.argv[2]
datadir = os.path.join(datadir, basedir)
resultsdir = exp.create_result_dir(os.path.join("multi-raxml", "norman", basedir + "_" + ranks))
result_msg = "multi-raxml git: \n" + exp.get_git_info(exp.multiraxml_root) + "\n"
result_msg += "raxml git: \n" + exp.get_git_info(exp.raxml_root) + "\n"
exp.write_results_info(resultsdir, result_msg) 
print("Results will be written in " + resultsdir) 


oldcommand_filename = os.path.join(datadir, "command.txt")
newcommand_filename = copyAndEditCommand(oldcommand_filename, resultsdir)
runCommand(newcommand_filename, resultsdir, ranks)


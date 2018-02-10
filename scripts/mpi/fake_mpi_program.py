import os
import sys
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp

if (len(sys.argv) != 4):
    print("Syntax error : ")
    print("python fake_mpi_program.py trees_number sites_number ranks")
    sys.exit(0)

trees = sys.argv[1]
sites = sys.argv[2]
ranks = sys.argv[3]
executable = os.path.join(exp.programs_root,
        "mpi",
        "fake_mpi_program",
        "build",
        "fake_mpi_program")

if (not os.path.isfile(executable)):
    print("Error: executable " + executable + " does not exist")
    sys.exit(0)

commands = []
commands.append("mpirun")
commands.append("-n")
commands.append(ranks)
commands.append(executable)
commands.append(trees)
commands.append(sites)
print("Running " + " ".join(commands))
subprocess.check_call(commands)

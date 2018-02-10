
import sys
import os
sys.path.insert(0, 'scripts')
import experiments as exp

# the goal of this script is to check that
# a haswell script that calls a python script
# that calls mpirun on several nodes
# will use all the nodes provided by haswell

# the answer is yes

if (len(sys.argv) != 4):
    print("Syntax error : ")
    print("python haswell_mpirun_from_python.py trees_number sites_number ranks")
    sys.exit(0)

trees = sys.argv[1]
sites = sys.argv[2]
ranks = sys.argv[3]


result_dir = exp.create_result_dir(os.path.join("mpi", "haswell_mpirun_from_python_" + ranks + "ranks"))
print("Results directory: " + result_dir)

submit_path = os.path.join(result_dir, "submit.sh")
script_path = os.path.join(exp.scripts_root, "mpi", "fake_mpi_program.py")
command = "python " + script_path + " " + trees + " " + sites + " " + ranks 

exp.submit_haswell(submit_path, command, ranks)

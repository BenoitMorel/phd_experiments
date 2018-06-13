import os
import experiments
import sys




if (len(sys.argv) != 5):
  print("Syntax: python scripts/launcher.py command_file cluster cores output_dir")
  print("cluster can be haswell or magny")
  sys.exit(1)


command_file = sys.argv[1]
cluster = sys.argv[2]
cores = int(sys.argv[3])
output_dir = sys.argv[4]

output_dir = os.path.join(experiments.results_root, "launcher", output_dir + "_" + str(cores), "run")
output_dir = experiments.create_result_dir(output_dir)

submit_file = os.path.join(output_dir, "launcher_submit.sh")
print("submit file" + submit_file)
command = "output_dir=" + output_dir + "\n" 
command += "cores=" + str(cores) + "\n"
command += "".join(open(command_file).readlines())

if (cluster == "haswell"):
  experiments.submit_haswell(submit_file, command, cores)
elif (cluster == "magny"):
  experiments.submit_magny(submit_file, command, cores)
else :
  print("unknown cluster value " + cluster)

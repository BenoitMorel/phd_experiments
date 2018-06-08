import os
import sys
import math

def submit_supermuc(submit_script, command, nodes):
  print("On " + str(nodes) + " nodes, from script " + submit_script + ", submit command " + command)
  
  directory = os.path.dirname(submit_script)
  nodes = int(nodes)
  with open(submit_script, "w") as f:
    f.write("#/bin/bash\n")
    f.write("#\n")
    f.write("#@ job_type = parallel\n")
    if (nodes < 20):
      f.write("#@ class = micro\n")
    else:
      f.write("#@ class = general\n")
    f.write("#\n")
    f.write("#@ node = " + str(nodes) + " \n")
    f.write("#\n")
    f.write("#@ tasks_per_node = 28\n")
    f.write("#@ wall_clock_limit = 24:00:00\n")
    f.write("#@ island_count = 1\n")
    f.write("#@ energy_policy_tag = examl3_aa_GAMMA\n")
    f.write("#@ minimize_time_to_solution = yes\n")
    f.write("## other example\n")
    f.write("#@ job_name = modeltest_ALL_1.1.aa_gamma\n")
    f.write("#@ network.MPI = sn_all,not_shared,us\n")
    f.write("#@ initialdir = " + directory + "\n")
    f.write("#@ output = job$(jobid).out\n")
    f.write("#@ error = job$(jobid).err\n")
    f.write("#@ notification=always\n")
    f.write("#@ notify_user=benoit.morel@h-its.org\n")
    f.write("#@ queue\n")
    f.write(". /etc/profile\n")
    f.write(". /etc/profile.d/modules.sh\n")
    f.write("\n")
    f.write(command)
  current_path = os.getcwd()
  os.chdir(directory)
  os.system("llsubmit " + os.path.basename(submit_script))
  os.chdir(current_path)


if (len(sys.argv) != 6):
  print("Syntax: python modeltest_split_and_run.py phyfile partfile splits per_splits_nodes outputdir")
  sys.exit(1)

# get arguments
phyfile = sys.argv[1]
partfile = sys.argv[2]
splits = int(sys.argv[3])
per_split_nodes = int(sys.argv[4])
outputdir = sys.argv[5]

# perform checks
if os.path.exists(outputdir):
  print("Directory " + outputdir + " already exists ")
  sys.exit(1)

# create output directory
os.makedirs(outputdir)
outputdir = os.path.abspath(outputdir)

print(outputdir)

modeltest = "/gpfs/work/pr58te/di68joy/1kite/software/benoit-modeltest/mpi/modeltest-mpi"
constant_arguments = " -r 1 -d aa -t mp " # seed = 1 type = aa tree = maximum parsimony


partitions = open(partfile, "r").read().split('\n')
partitions_count = len(partitions)
split_len = math.ceil(partitions_count / splits)
for i in range(0, splits):
  split_outputdir = os.path.join(outputdir, "split_" + str(i))
  os.makedirs(split_outputdir)
  split_partition = os.path.join(split_outputdir, "split_partition.part")
  with open(split_partition, "w") as f:
    begin = i * split_len
    end = (i + 1) * split_len
    f.write('\n'.join(partitions[begin:end]))
  split_script = os.path.join(split_outputdir, "MT_split_" + str(i) + ".sh")
  command = "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/gpfs/work/pr58te/di68joy/1kite/software/benoit-modeltest/build/lib\n"
  command += "poe " + modeltest + constant_arguments  + " -q " + split_partition
  command += " -i " + phyfile
  command += " -o " + split_outputdir 
  submit_supermuc(split_script, command, per_split_nodes)



import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp

def generate_scheduler_commands_file(dataset_dir, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = os.path.join(dataset_dir, "speciesTree.newick")
  
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      treerecs_dir = os.path.join(family_dir, "treerecs")
      try:
        os.makedirs(treerecs_dir)
      except:
        pass
      
      alignment_descriptor = os.path.join(treerecs_dir, "alignment_descriptor.txt")
      with open(alignment_descriptor, "w") as ali_writer:
        ali_writer.write("GTR\n" + os.path.abspath(os.path.join(family_dir, "alignment.msa")))
      treerecs_output = os.path.join(treerecs_dir, "treerecs_output")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(exp.treerecs_exec)
      command.append("--seed")
      command.append("42")
      command.append("-g")
      command.append(os.path.join(family_dir, "raxmlGeneTree.newick"))
      command.append("-s")
      command.append(speciesTree)
      command.append("-o")
      command.append(treerecs_output)
      command.append("-a")
      command.append(alignment_descriptor)
      command.append("-t")
      command.append("all")
      command.append("--ale-evaluation")
      command.append("-T")
      command.append("7")
      command.append("--select-best-tree")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def generate_scheduler_command(command_file, cores, output_dir):
  command = ""
  parallelization = "onecore"
  command += "mpirun -np " + str(cores) + " "
  command += exp.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += exp.treerecs_exec + " "
  command += command_file + " "
  command += output_dir + " " 
  command += "0"
  return command 

def extract_treerecs_trees(families_dir):
  for family in os.listdir(families_dir):
    treerecsTree = os.path.join(families_dir, family, "treerecs", "treerecs_output.newick.best")
    if (os.path.isfile(treerecsTree)):
      lines = open(treerecsTree).readlines()
      with open(os.path.join(families_dir, family, "treerecsGeneTree.newick"), "w") as writer:
        for line in lines:
          if (not line.startswith(">")):
            writer.write(line)
    else:
      print("Warning: no treerecs tree for family " + family)
  

def run_treerecs_on_families(dataset_dir, cores):
  output_dir = os.path.join(dataset_dir, "treerecs_run")
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(dataset_dir, cores, output_dir)
  command = generate_scheduler_command(scheduler_commands_file, cores, output_dir)
  print(command.split(" "))
  subprocess.check_call(command.split(" "))
  extract_treerecs_trees(os.path.join(dataset_dir, "families"))
  
if (__name__== "__main__"):
  max_args_number = 3
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_treerecs dataset_dir cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  cores = int(sys.argv[2])
  run_treerecs_on_families(dataset_dir, cores)



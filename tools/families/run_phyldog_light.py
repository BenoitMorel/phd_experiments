import os
import sys
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp


def convertToPhyldogSpeciesTree(speciesTree, phyldogSpeciesTree):
  command = "sed s/)n[0123456789]*/)/g " + speciesTree #+ " > " + phyldogSpeciesTree
  print(command.split(" "))
  with open(phyldogSpeciesTree, "w") as output:
    subprocess.check_call(command.split(" "), stdout=output)

def generate_scheduler_commands_file(dataset_dir, cores, output_dir, scheduler_output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(scheduler_output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = os.path.join(dataset_dir, "speciesTree.newick")
  phyldogSpeciesTree = os.path.join(dataset_dir, "phyldogSpeciesTree.newick")
  convertToPhyldogSpeciesTree(speciesTree, phyldogSpeciesTree)

  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      phyldog_dir = os.path.join(family_dir, "phyldog")
      try:
        os.makedirs(phyldog_dir)
      except:
        pass
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("species.tree.file=" + phyldogSpeciesTree)
      command.append("gene.tree.file=" + os.path.join(family_dir, "raxmlGeneTree.newick"))
      command.append("input.sequence.file=" + os.path.join(family_dir, "alignment.msa"))
      command.append("taxaseq.file=" + os.path.join(family_dir, "mapping.link"))
      command.append("likelihood.evaluator=LIBPLL2")
      command.append("model=GTR ")
      os.makedirs(os.path.join(results_dir, family))
      command.append("output.file=" + os.path.join(family_dir, "phyldog", "phyldog"))
      #command.append("output.file=" + os.path.join(results_dir, family, "phyldog"))
      #print(" ".join(command))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def generate_scheduler_command(command_file, cores, scheduler_output_dir):
  command = ""
  parallelization = "onecore"
  command += "mpirun -np " + str(cores) + " "
  command += exp.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += exp.phyldog_light_exec + " "
  command += command_file + " "
  command += scheduler_output_dir + " " 
  command += "0"
  return command 

def run_phyldog_light_on_families(dataset_dir, cores):
  output_dir = os.path.join(dataset_dir, "phyldog_run")
  os.makedirs(output_dir)
  scheduler_output_dir = os.path.join(output_dir, "scheduler")
  scheduler_commands_file = generate_scheduler_commands_file(dataset_dir, cores, output_dir, scheduler_output_dir)
  command = generate_scheduler_command(scheduler_commands_file, cores, scheduler_output_dir)
  print(command.split(" "))
  subprocess.check_call(command.split(" "))


if (__name__== "__main__"):
  max_args_number = 3
  if len(sys.argv) < max_args_number:
    print("Syntax error: python families_phyldog_light_launcher.py dataset_dir cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  cores = int(sys.argv[2])
  run_phyldog_light_on_families(dataset_dir, cores)


import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "notung"))
import experiments as exp
import events_scenario_extraction as extract
import convert_to_notung_tree
  
def generate_notung_files(dataset_dir):
  families_dir = os.path.join(dataset_dir, "families")
  for family in os.listdir(families_dir):
    family_dir = os.path.join(families_dir, family)
    input_tree = os.path.join(family_dir, "raxmlGeneTree.newick")
    input_species_tree = os.path.join(family_dir, "speciesTree.newick")
    mapping_file = os.path.join(family_dir, "mapping.link")
    notung_tree = os.path.join(family_dir, "raxmlGeneTree.notung")
    notung_species_tree = os.path.join(family_dir, "speciesTree.notung")
    notung_mapping = os.path.join(family_dir, "mapping.notung")
    convert_to_notung_tree.convert_to_notung_tree(input_tree, input_species_tree, mapping_file, notung_tree, notung_species_tree, notung_mapping)

def back_convert_notung_files(dataset_dir):
  families_dir = os.path.join(dataset_dir, "families")
  for family in os.listdir(families_dir):
    family_dir = os.path.join(families_dir, family)
    notung_tree = os.path.join(family_dir, "raxmlGeneTree.notung.rearrange.0")
    notung_mapping = os.path.join(family_dir, "mapping.notung")
    output_tree = os.path.join(family_dir, "notungGeneTree.newick")
    convert_to_notung_tree.back_convert_notung_tree(notung_tree, notung_mapping, output_tree)

def generate_scheduler_command(command_file, cores, output_dir):
  command = ""
  parallelization = "onecore"
  command += "mpirun -np " + str(cores) + " "
  command += exp.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += "java" + " "
  command += command_file + " "
  command += output_dir + " " 
  command += "0"
  return command 
  
def generate_scheduler_commands_file(dataset_dir, threshold, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = os.path.join(dataset_dir, "speciesTree.notung")
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      #command.append("java")
      command.append("-jar")
      command.append(exp.notung_jar)
      command.append(os.path.join(family_dir, "raxmlGeneTree.notung"))
      command.append("-s")
      command.append(os.path.join(family_dir, "speciesTree.notung"))
      command.append("--speciestag")
      command.append("postfix")
      command.append("--rearrange")
      command.append("--outputdir")
      command.append(family_dir)
      command.append("--edgeweights")
      command.append("name")
      command.append("--threshold")
      command.append(str(threshold))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def run_notung_on_families(dataset_dir, threshold, cores):
  output_dir = os.path.join(dataset_dir, "notung_run")
  os.makedirs(output_dir)
  generate_notung_files(dataset_dir)
  scheduler_commands_file = generate_scheduler_commands_file(dataset_dir,  threshold, cores, output_dir)
  command = generate_scheduler_command(scheduler_commands_file, cores, output_dir)
  print(command.split(" "))
  subprocess.check_call(command.split(" "))
  back_convert_notung_files(dataset_dir)


if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_notung dataset_dir threshold cores.")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  threshold = int(sys.argv[2])
  cores = int(sys.argv[3])

  run_notung_on_families(dataset_dir, threshold, cores)




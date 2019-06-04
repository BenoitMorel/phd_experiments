import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "notung"))
import experiments as exp
import convert_to_notung_tree
import time
import saved_metrics
  
def generate_notung_files(datadir):
  for family in fam.get_families_list(datadir):
    misc_dir = fam.get_misc_dir(datadir, family)
    input_tree = fam.get_raxml_tree(datadir, family) 
    input_species_tree = fam.get_species_tree(datadir)
    mapping_file = fam.get_mappings(datadir, family)
    notung_tree = os.path.join(misc_dir, "raxmlGeneTree.notung")
    notung_species_tree = os.path.join(misc_dir, "speciesTree.notung")
    notung_mapping = os.path.join(misc_dir, "mapping.notung")
    convert_to_notung_tree.convert_to_notung_tree(input_tree, input_species_tree, mapping_file, notung_tree, notung_species_tree, notung_mapping)

def back_convert_notung_files(datadir, threshold):
  for family in fam.get_families_list(datadir):
    misc_dir = fam.get_misc_dir(datadir, family)
    notung_tree = os.path.join(misc_dir, "raxmlGeneTree.notung.rearrange.0")
    notung_mapping = os.path.join(misc_dir, "mapping.notung")
    output_tree = fam.get_notung_tree(datadir, family, threshold)
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
  
def generate_scheduler_commands_file(datadir, threshold, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = os.path.join(datadir, "misc", "speciesTree.notung")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      misc_dir = fam.get_misc_dir(datadir, family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      #command.append("java")
      command.append("-jar")
      command.append(exp.notung_jar)
      command.append(os.path.join(misc_dir, "raxmlGeneTree.notung"))
      command.append("-s")
      command.append(os.path.join(misc_dir, "speciesTree.notung"))
      command.append("--speciestag")
      command.append("postfix")
      command.append("--rearrange")
      command.append("--outputdir")
      command.append(misc_dir)
      command.append("--edgeweights")
      command.append("name")
      command.append("--threshold")
      command.append(str(threshold))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file

def run_notung_on_families(datadir, threshold, cores):
  output_dir = os.path.join(datadir, "runs", "notung" + str(threshold) + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  generate_notung_files(datadir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir,  threshold, cores, output_dir)
  
  start = time.time()
  exp.run_with_scheduler("java", scheduler_commands_file, "onecore", cores, output_dir, "logs.txt")   
  saved_metrics.save_metrics(datadir, "notung" + str(int(threshold)), (time.time() - start), "runtimes") 
  back_convert_notung_files(datadir, threshold)


if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_notung datadir threshold cores.")
    sys.exit(0)


  datadir = sys.argv[1]
  threshold = int(sys.argv[2])
  cores = int(sys.argv[3])

  run_notung_on_families(datadir, threshold, cores)





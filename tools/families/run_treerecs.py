import os
import sys
import subprocess
import shutil
import time
import saved_metrics
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import fam

def generate_scheduler_commands_file(datadir, is_dna, cores, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = fam.get_species_tree(datadir)
  model = "GTR"
  if (not is_dna):
    model = "LG"
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      treerecs_dir = fam.get_misc_dir(datadir, family)
      alignment_descriptor = os.path.join(treerecs_dir, "alignment_descriptor.txt")
      with open(alignment_descriptor, "w") as ali_writer:
        ali_writer.write(model + "\n" + os.path.abspath(os.path.join(family_dir, "alignment.msa")))
      treerecs_output = os.path.join(treerecs_dir, "treerecs_output")
      mapping_file = fam.get_treerecs_mappings(datadir, family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append(exp.treerecs_exec)
      command.append("--seed")
      command.append("42")
      command.append("-g")
      command.append(fam.get_raxml_tree(datadir, family))
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
      if(os.path.isfile(mapping_file)):
        command.append("-S")
        command.append(mapping_file)

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
  print(command)
  return command 

def extract_treerecs_trees(datadir):
  families_dir = os.path.join(datadir, "families")
  for family in os.listdir(families_dir):
    treerecsTree = os.path.join(families_dir, family, "misc", "treerecs_output.newick.best")
    if (os.path.isfile(treerecsTree)):
      lines = open(treerecsTree).readlines()
      with open(fam.get_treerecs_tree(datadir, family), "w") as writer:
        for line in lines:
          if (not line.startswith(">")):
            writer.write(line)
    else:
      print("Warning: no treerecs tree for family " + family)  

def run_treerecs_on_families(datadir, is_dna, cores):
  output_dir = os.path.join(datadir, "runs", "treerecs_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, is_dna, cores, output_dir)
  command = generate_scheduler_command(scheduler_commands_file, cores, output_dir)
  print(command.split(" "))
  start = time.time()
  subprocess.check_call(command.split(" "), stdout = sys.stdout)
  saved_metrics.save_metrics(datadir, "Treerecs", (time.time() - start), "runtimes") 
  extract_treerecs_trees(datadir)
  
if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_treerecs datadir is_dna cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  datadir = sys.argv[1]
  is_dna = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])

  run_treerecs_on_families(datadir, is_dna, cores)


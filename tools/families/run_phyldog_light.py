import os
import sys
import subprocess
import shutil
import time
import saved_metrics
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import events_scenario_extraction as extract
import experiments as exp


def convertToPhyldogSpeciesTree(speciesTree, phyldogSpeciesTree):
  command = "sed s/)[nH][0123456789]*/)/g " + speciesTree #+ " > " + phyldogSpeciesTree
  print(command.split(" "))
  with open(phyldogSpeciesTree, "w") as output:
    subprocess.check_call(command.split(" "), stdout=output)

def generate_scheduler_commands_file(dataset_dir, is_dna, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = os.path.join(dataset_dir, "speciesTree.newick")
  phyldogSpeciesTree = os.path.join(dataset_dir, "phyldogSpeciesTree.newick")
  convertToPhyldogSpeciesTree(speciesTree, phyldogSpeciesTree)
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      phyldog_dir = os.path.join(family_dir, "misc")
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("species.tree.file=" + phyldogSpeciesTree)
      command.append("gene.tree.file=" + fam.getRaxmlTree(dataset_dir, family))
      command.append("input.sequence.file=" + os.path.join(family_dir, "alignment.msa"))
      command.append("taxaseq.file=" + os.path.join(family_dir, "mapping.link"))
      if (is_dna):
        command.append("model=GTR")
        command.append("likelihood.evaluator=LIBPLL2")
      else:
        command.append("model=LG08")
        command.append("likelihood.evaluator=PLL")
      if (not is_dna):
        command.append("alphabet=Protein")
      os.makedirs(os.path.join(results_dir, family))
      command.append("output.file=" + os.path.join(phyldog_dir, "phyldog"))
      command.append("branch.expected.numbers.optimization=average")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def generate_scheduler_command(command_file, cores, output_dir):
  command = ""
  parallelization = "onecore"
  command += "mpirun -np " + str(cores) + " "
  command += exp.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += exp.phyldog_light_exec + " "
  command += command_file + " "
  command += output_dir + " " 
  command += "0"
  return command 

def extract_phyldog_trees(dataset_dir):
  families_dir = os.path.join(dataset_dir, "families")
  for family in os.listdir(families_dir):
    phyldogTree = os.path.join(families_dir, family, "misc", "phyldog_reconciled.tree")
    if (os.path.isfile(phyldogTree)):
      shutil.copyfile(phyldogTree, fam.getPhyldogTree(dataset_dir, family))
    else:
      print("Warning: no phyldog tree for family " + family)
  

def run_phyldog_light_on_families(dataset_dir, is_dna, cores):
  fam.init_dataset_dir(dataset_dir)
  output_dir = os.path.join(dataset_dir, "runs", "phyldog_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(dataset_dir, is_dna, cores, output_dir)
  command = generate_scheduler_command(scheduler_commands_file, cores, output_dir)
  print(command.split(" "))
  start = time.time()
  subprocess.check_call(command.split(" "), stdout = sys.stdout)
  saved_metrics.save_metrics(dataset_dir, "Phyldog", (time.time() - start), "runtimes") 
  extract_phyldog_trees(dataset_dir)
  extract.extract_events_from_phyldog(dataset_dir)

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_phyldog_light.py dataset_dir is_dna cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  is_dna = int(sys.argv[2]) != 0
  cores = int(sys.argv[3])
  run_phyldog_light_on_families(dataset_dir, is_dna, cores)


import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
import experiments as exp


def mymakedirs(path):
  try:
    os.makedirs(path)
  except:
    pass

def convert_to_phyldog_species_tree(speciesTree, phyldogSpeciesTree):
  command = "sed s/)n[0123456789]*/)/g " + speciesTree #+ " > " + phyldogSpeciesTree
  print(command.split(" "))
  with open(phyldogSpeciesTree, "w") as output:
    subprocess.check_call(command.split(" "), stdout=output)

def generate_scheduler_commands_file(dataset_dir, cores, output_dir):
  families_dir = os.path.join(dataset_dir, "families")
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  speciesTree = os.path.join(dataset_dir, "speciesTree.newick")
  phyldogSpeciesTree = os.path.join(dataset_dir, "phyldogSpeciesTree.newick")
  convert_to_phyldog_species_tree(speciesTree, phyldogSpeciesTree)
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_dir = os.path.join(families_dir, family)
      phyldog_dir = os.path.join(family_dir, "phyldog")
      mymakedirs(phyldog_dir)
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
      mymakedirs(os.path.join(results_dir, family))
      command.append("output.file=" + os.path.join(family_dir, "phyldog", "phyldog"))
      #command.append("output.file=" + os.path.join(results_dir, family, "phyldog"))
      #print(" ".join(command))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file
     
def extract_phyldog_trees(families_dir):
  for family in os.listdir(families_dir):
    phyldogTree = os.path.join(families_dir, family, "phyldog", "phyldog_reconciled.tree")
    if (os.path.isfile(phyldogTree)):
      shutil.copyfile(phyldogTree, os.path.join(families_dir, family, "phyldogGeneTree.newick"))
    else:
      print("Warning: no phyldog tree for family " + family)
  

def run_phyldog_on_families(dataset_dir, cores):
  output_dir = os.path.join(dataset_dir, "phyldog_run")
  mymakedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(dataset_dir, cores, output_dir)
  exp.run_with_scheduler(scheduler_commands_file, "onecore", cores, output_dir)   
  extract_phyldog_trees(os.path.join(dataset_dir, "families"))

if (__name__== "__main__"):
  max_args_number = 3
  if len(sys.argv) < max_args_number:
    print("Syntax error: python run_phyldog.py dataset_dir cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)


  dataset_dir = sys.argv[1]
  cores = int(sys.argv[2])
  run_phyldog_on_families(dataset_dir, cores)


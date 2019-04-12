import os
import sys
import subprocess
import shutil
import time
import saved_metrics
import fam
import time
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
      command.append("likelihood.evaluator=PLL")
      if (is_dna):
        command.append("model=GTR")
      else:
        command.append("model=LG08")
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

def clean_phyldog(dataset_dir):
  dataset_name = os.path.basename(dataset_dir)
  print("CLEANING " + dataset_name)
  for f in os.listdir("."):
    if (f.startswith("tmpPLL") and  (dataset_name.replace(".", "_") in f.replace(".", "_"))):
      os.remove(f)

def get_phyldog_run_dir(dataset_dir):
  return os.path.join(dataset_dir, "runs", "phyldog_run")

def run_phyldog_light_on_families(dataset_dir, is_dna, cores):
  fam.init_dataset_dir(dataset_dir)
  output_dir = get_phyldog_run_dir(dataset_dir)
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  generate_options(dataset_dir, is_dna)
  run_phyldog(dataset_dir, cores)
  extract_phyldog(dataset_dir)
  return
  scheduler_commands_file = generate_scheduler_commands_file(dataset_dir, is_dna, cores, output_dir)
  command = generate_scheduler_command(scheduler_commands_file, cores, output_dir)
  print(command.split(" "))
  start = time.time()
  subprocess.check_call(command.split(" "), stdout = sys.stdout)
  saved_metrics.save_metrics(dataset_dir, "Phyldog", (time.time() - start), "runtimes") 
  clean_phyldog(dataset_dir)
  extract_phyldog_trees(dataset_dir)
  extract.extract_events_from_phyldog(dataset_dir)
    
    
def add_starting_tree(option_file, tree_path):
  lines = open(option_file).readlines()
  with open(option_file, "w") as writer:
    for line in lines:
      if (line.startswith("init.gene.tree")):
        continue
      if (line.startswith("gene.tree.file")):
          continue
      writer.write(line)
    writer.write("\n")
    writer.write("init.gene.tree=user\n")
    writer.write("gene.tree.file=" + tree_path + "\n")

def generate_options(dataset_dir, is_dna):
  families_dir = os.path.join(dataset_dir, "families")
  speciesTree = os.path.join(dataset_dir, "speciesTree.newick")
  phyldogSpeciesTree = os.path.join(dataset_dir, "phyldogSpeciesTree.newick")
  convertToPhyldogSpeciesTree(speciesTree, phyldogSpeciesTree)
  phyldog_run_dir = get_phyldog_run_dir(dataset_dir)
  prepare_input = os.path.join(phyldog_run_dir, "prepare_input.txt")
  dataset_dir = os.path.abspath(dataset_dir)
  phyldog_run_dir = os.path.abspath(phyldog_run_dir)
  mappings_dir = os.path.join(phyldog_run_dir, "mappings")
  exp.reset_dir(mappings_dir)
  all_raxml_trees_dir = os.path.join(phyldog_run_dir, "all_raxml_trees")
  exp.reset_dir(all_raxml_trees_dir)
  for family in os.listdir(families_dir):
    family_dir = os.path.join(families_dir, family)
    phyldog_mapping = os.path.join(family_dir, "mapping.link")
    new_phyldog_mapping = os.path.join(mappings_dir, family + ".link")
    exp.relative_symlink(phyldog_mapping, new_phyldog_mapping)
    old_raxml_tree = fam.getRaxmlTree(dataset_dir, family)
    new_raxml_tree = os.path.join(all_raxml_trees_dir, family + ".newick")
    exp.relative_symlink(old_raxml_tree, new_raxml_tree) 
  with open(prepare_input, "w") as writer:
    writer.write(os.path.join(dataset_dir, "alignments") + "\n")
    if (is_dna):
      writer.write("DNA\n" )
    else:
      writer.write("PROTEIN\n")
    writer.write("FASTA\n")
    writer.write(os.path.join(phyldog_run_dir, "mappings") + "\n")
    writer.write(os.path.join(phyldog_run_dir, "options") + "\n")
    writer.write(os.path.join(phyldog_run_dir, "results") + "\n")
    writer.write("yes" + "\n")
    writer.write(os.path.join(dataset_dir, "phyldogSpeciesTree.newick") + "\n")
    writer.write("no" + "\n") #opt species tree
    writer.write("yes" + "\n") # opt dup loss
    writer.write("branchwise" + "\n") # branchwise DL opt
    writer.write("no" + "\n") # same number of genes ?
    writer.write("yes" + "\n") # opt gene trees
    writer.write("48" + "\n") # max time (hours)
  prepare_data_script = os.path.join(exp.tools_root, "families", "prepareData.py")
  subprocess.check_call(["/bin/bash", "-c", "python " + prepare_data_script + " < " + prepare_input])
  time.sleep(1)
  for family in os.listdir(families_dir):
    option_file = os.path.join(phyldog_run_dir, "options", family + ".opt")
    add_starting_tree(option_file, fam.getRaxmlTree(dataset_dir, family))


   
def run_phyldog(dataset_dir, cores):
  families_number = len(os.listdir(os.path.join(dataset_dir, "families")))
  print("plop")
  print(str(families_number) + " families")
  cores = str(min(families_number, int(cores)))
  phyldog_run_dir = get_phyldog_run_dir(dataset_dir)
  command = []
  command.append("mpirun")
  command.append("-n")
  command.append(str(cores))
  command.append(exp.phyldog_exec)
  command.append("param=" + os.path.join(phyldog_run_dir, "options", "GeneralOptions.txt"))
  print("EXECUTE PHYLDOG")
  print(" ".join(command))
  logs = open(os.path.join(phyldog_run_dir, "logs.txt"), "w")
  subprocess.check_call(command, stdout = logs)
  print("end EXECUTE PHYLDOG")


def extract_phyldog(dataset_dir):
  results_dir = os.path.join(get_phyldog_run_dir(dataset_dir), "results")
  families_dir = os.path.join(dataset_dir, "families")
  for family in os.listdir(families_dir):
    phyldog_tree = os.path.join(results_dir, family + ".ReconciledTree")
    new_phyldog_tree = fam.getPhyldogTree(dataset_dir, family)
    shutil.copy(phyldog_tree, new_phyldog_tree)


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


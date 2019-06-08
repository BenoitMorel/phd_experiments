import sys
import os
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/msa_edition')
sys.path.insert(0, 'tools/families')
import analyze_dataset
import experiments as exp
import exp_jointsearch_utils as utils
import msa_analyzer
import subprocess
import shutil
import fam

parallelization = "split"
nodes_per_core = 2

def get_generax_datasets():
  root_datadir = os.path.join(exp.benoit_datasets_root, "families")
  datasets = {}
  for dataset in os.listdir(root_datadir):
    datasets[dataset] = os.path.join(root_datadir, dataset)
  return datasets
datasets = get_generax_datasets()

def generate_scheduler_commands_file(datadir, starting_tree, strategy, nodes_per_core, additional_arguments, cores, output_dir, scheduler_output_dir):
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  results_dir = os.path.join(scheduler_output_dir, "results")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_path = fam.get_family_path(datadir, family)
      gene_tree = fam.get_gene_tree(family_path, starting_tree)
      if (gene_tree != "__random__" and len(open(gene_tree).readlines()) == 0):
        continue
      alignment =  fam.get_alignment(datadir, family)
      tree_size = msa_analyzer.get_taxa_number(alignment)
      if (0 == nodes_per_core):
        writer.write(family + " 1 " + str(tree_size) + " ")
      else:
        max_cores = 16
        cores = max(min(int(tree_size / nodes_per_core), max_cores), 1)
        writer.write(family + " " + str(cores) + " " + str(tree_size) + " ")
      writer.write("-g " + gene_tree + " ")
      writer.write("-s " + fam.get_species_tree(datadir) +  " ")
      writer.write("-m " + fam.get_mappings(datadir, family) + " ")
      writer.write("-a " + alignment +  " ")
      os.makedirs(os.path.join(results_dir, family))
      writer.write("-p " + os.path.join(results_dir, family, "jointsearch") + " ") 
      writer.write("--strategy " + strategy + " ")
      writer.write(" ".join(additional_arguments) + " ")
      writer.write("\n")

  return scheduler_commands_file


def extract_trees(datadir, results_family_dir, prefix):
  print("run name " + prefix)
  results_dir = os.path.join(results_family_dir, "results")
  for family in fam.get_families_list(datadir):
    family_results = fam.get_results(datadir, family)
    source = os.path.join(results_dir, family, "jointsearch.newick")
    dest = os.path.join(family_results, prefix + "GeneTree.newick")
    shutil.copy(source, dest)


def launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments):
  command = ["python"]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("MultipleJointSearch", dataset, strategy + "_start_" + starting_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  result_msg = "GeneRax git: \n" + exp.get_git_info(exp.joint_search_root)
  exp.write_results_info(resultsdir, result_msg) 
  submit_path = os.path.join(resultsdir, "submit.sh")
  command.append(resultsdir)
  exp.submit(submit_path, " ".join(command), cores, cluster) 

def get_exec(parallelization):
  if (parallelization == "split"):
    return exp.joint_search_lib
  elif (parallelization == "openmp" or parallelization == "onecore" or parallelization == "fork"):
    return exp.joint_search_exec 
  else:
    print("Unknown scheduler parallelization " + parallelization)
    sys.exit(0)

def run(dataset, strategy, starting_tree, cores, additional_arguments, resultsdir):
  is_protein = exp.checkAndDelete("--protein", additional_arguments)
  run_name = exp.getAndDelete("--run", additional_arguments, "lastRun") 
  #mode = get_mode_from_additional_arguments(additional_arguments)
  datadir = datasets[dataset]
  families_dir = os.path.join(datadir, "families")
  scheduler_output_dir = os.path.join(resultsdir, "scheduler_run")
  os.makedirs(scheduler_output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, starting_tree, strategy, nodes_per_core, additional_arguments, cores, resultsdir, scheduler_output_dir)
  exp.run_with_scheduler(get_exec(parallelization), scheduler_commands_file, parallelization, cores,  scheduler_output_dir)
  extract_trees(datadir, os.path.join(resultsdir, "scheduler_run"), run_name)
  analyze_dataset.analyze(datadir, run_name)
  print("Output in " + resultsdir)

if (__name__ == "__main__"): 
  is_run = ("--exprun" in sys.argv)
  resultsdir = ""
  if (is_run):
    resultsdir = sys.argv[-1]
    sys.argv = sys.argv[:-2]
  
  min_args_number = 6
  min_args_number = 6
  if (len(sys.argv) < min_args_number):
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset strategy starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
    for dataset in datasets:
      print("\t" + dataset)
    sys.exit(1)
  
  dataset = sys.argv[1]
  strategy = sys.argv[2]
  starting_tree = sys.argv[3]
  cluster = sys.argv[4]
  cores = int(sys.argv[5])
  additional_arguments = sys.argv[min_args_number:]

  if (is_run):
    run(dataset, strategy, starting_tree, cores, additional_arguments, resultsdir)
  else:
    launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments)



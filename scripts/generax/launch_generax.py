import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import saved_metrics
import analyze_dataset
import experiments as exp
import shutil
import time
import fam

def get_possible_strategies():
  return ["SPR", "EVAL"]

def check_inputs(strategy):
  if (not (strategy in get_possible_strategies())):
    print("Unknown search strategy " + strategy)
    exit(1)


def get_generax_datasets():
  root_datadir = os.path.join(exp.benoit_datasets_root, "families")
  datasets = {}
  for dataset in os.listdir(root_datadir):
    datasets[dataset] = os.path.join(root_datadir, dataset)
  return datasets

datasets = get_generax_datasets()

def build_generax_families_file(dataset, starting_tree, is_protein, output):
  families_dir = os.path.join(dataset, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    plop = 0
    for family in os.listdir(families_dir):
      
      family_path = os.path.join(families_dir, family)
      writer.write("- " + family + "\n")
      writer.write("starting_gene_tree = " + fam.get_gene_tree(family_path, starting_tree) + "\n")
      writer.write("alignment = " + fam.get_alignment_file(family_path) + "\n")
      writer.write("mapping = " + fam.get_mappings(dataset, family) + "\n")
      raxml_model = ""
      if (starting_tree != "random"):
        raxml_model = fam.get_raxml_model(family_path)
      if (os.path.isfile(raxml_model)):
        writer.write("subst_model = " + raxml_model + "\n")
      else:
        if (is_protein):
          writer.write("subst_model = LG\n")
        else:
          writer.write("subst_model = GTR\n")

def get_generax_command(generax_families_file, species_tree, strategy, additional_arguments, output_dir, mode, cores):
    executable = exp.generax_exec
    if (mode == "gprof"):
      executable = exp.generax_gprof_exec
    elif (mode == "scalasca"):
      executable = exp.generax_scalasca_exec
    generax_output = os.path.join(output_dir, "generax")
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-f")
    command.append(generax_families_file)
    command.append("-s")
    command.append(species_tree)
    command.append("--strategy")
    command.append(strategy)
    command.append("-p")
    command.append(generax_output)
    command.extend(additional_arguments)
    return " ".join(command)

def run_generax(datadir, strategy, generax_families_file, mode, cores, additional_arguments, resultsdir):
  species_tree = fam.get_species_tree(datadir)
  command = get_generax_command(generax_families_file, species_tree, strategy, additional_arguments, resultsdir, mode, cores)
  print(command)
  subprocess.check_call(command.split(" "), stdout = sys.stdout)


def get_mode_from_additional_arguments(additional_arguments):
  mode = "normal"
  if ("--scalasca" in additional_arguments):
    mode = "scalasca"
    additional_arguments.remove("--scalasca")
  elif ("--gprof" in additional_arguments):
    mode = "gprof"
    additional_arguments.remove("--gprof")
  return mode


def extract_trees(data_family_dir, results_family_dir, prefix):
  results_dir = os.path.join(results_family_dir, "results")
  for msa in os.listdir(results_dir):
    output_msa_dir = os.path.join(data_family_dir, msa, "results")
    try:
      os.mkdir(output_msa_dir)
    except:
      pass
    source = os.path.join(results_dir, msa, "geneTree.newick")
    dest = os.path.join(output_msa_dir, prefix + "GeneTree.newick")
    shutil.copy(source, dest)


def run(dataset, strategy, starting_tree, cores, additional_arguments, resultsdir):
  is_protein = exp.checkAndDelete("--protein", additional_arguments)
  run_name = exp.getAndDelete("--run", additional_arguments, "lastRun") 
  mode = get_mode_from_additional_arguments(additional_arguments)
  if (not dataset in datasets):
    print("Error: " + dataset + " is not in " + str(datasets))
    exit(1)
  datadir = datasets[dataset]
  generax_families_file = os.path.join(resultsdir, "generax_families.txt")
  build_generax_families_file(datadir, starting_tree, is_protein, generax_families_file)
  start = time.time()
  run_generax(datadir, strategy, generax_families_file, mode, cores, additional_arguments, resultsdir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
  extract_trees(os.path.join(datadir, "families"), os.path.join(resultsdir, "generax"), run_name)
  try:
    analyze_dataset.analyze(datadir, run_name)
  except:
    print("Analyze failed!!!!")

  print("Output in " + resultsdir)

def launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments):
  command = ["python"]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("GeneRax", dataset, strategy + "_start_" + starting_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  #result_msg = "GeneRax git: \n" + exp.get_git_info(exp.joint_search_root)
  #result_msg = ""
  #exp.write_results_info(resultsdir, result_msg) 
  submit_path = os.path.join(resultsdir, "submit.sh")
  command.append(resultsdir)
  exp.submit(submit_path, " ".join(command), cores, cluster) 
  

if (__name__ == "__main__"): 
  is_run = ("--exprun" in sys.argv)
  resultsdir = ""
  if (is_run):
    resultsdir = sys.argv[-1]
    sys.argv = sys.argv[:-2]
    
  min_args_number = 6
  if (len(sys.argv) < min_args_number):
    for dataset in datasets:
      print("\t" + dataset)
    print("strategy: " + ",".join(get_possible_strategies()))
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset strategy starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
    sys.exit(1)

  dataset = sys.argv[1]
  strategy = sys.argv[2]
  starting_tree = sys.argv[3]
  cluster = sys.argv[4]
  cores = int(sys.argv[5])
  additional_arguments = sys.argv[min_args_number:]
  check_inputs(strategy)

  if (starting_tree == "raxml"):
    print("use raxml-ng instead of raxml please")
    exit(1)

  if (is_run):
    run(dataset, strategy, starting_tree, cores, additional_arguments, resultsdir)
  else:
    launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments)




import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import saved_metrics
import rf_cells
import experiments as exp
import shutil
import time
import fam
import sequence_model
import fast_rf_cells

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


def has_multiple_sample(starting_tree):
  return "ale" in starting_tree.lower() or "multiple" in starting_tree.lower()

def get_starting_tree_path(datadir, subst_model, family, starting_tree):
  if (has_multiple_sample(starting_tree)):
    return os.path.join(fam.get_family_misc_dir(datadir, family), starting_tree + "." + subst_model + "_onesample.geneTree")
  else:
    return fam.build_gene_tree_path(datadir, subst_model, family, starting_tree)

# GeneRax does not accept tree files with multiple trees
def sample_one_starting_tree(datadir, subst_model, starting_tree):
  for family in fam.get_families_list(datadir):
    input_tree = fam.build_gene_tree_path(datadir, subst_model, family, starting_tree)
    output_tree = get_starting_tree_path(datadir, subst_model, family, starting_tree)
    tree = open(input_tree, "r").readline()
    open(output_tree, "w").write(tree)

def build_generax_families_file(datadir, starting_tree, subst_model, output):
  if (has_multiple_sample(starting_tree)):
    sample_one_starting_tree(datadir, subst_model, starting_tree)
  families_dir = os.path.join(datadir, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    plop = 0
    print("starting gene tree " + starting_tree)
    for family in os.listdir(families_dir):
      
      family_path = os.path.join(families_dir, family)
      writer.write("- " + family + "\n")
      gene_tree = get_starting_tree_path(datadir, subst_model, family, starting_tree)
      if (starting_tree == "random"):
        gene_tree = "__random__"
      writer.write("starting_gene_tree = " + gene_tree + "\n")
      writer.write("alignment = " + fam.get_alignment_file(family_path) + "\n")
      writer.write("mapping = " + fam.get_mappings(datadir, family) + "\n")
      raxml_model = ""
      if (starting_tree != "random" and starting_tree != "true"):
        raxml_model = fam.get_raxml_best_model(datadir, subst_model, family)
      if (os.path.isfile(raxml_model)):
        writer.write("subst_model = " + raxml_model + "\n")
      else:
        writer.write("subst_model = " + sequence_model.get_raxml_model(subst_model) + "\n")

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
  command = ["mpirun", "-np", str(cores), "sleep", "5"]
  subprocess.check_call(command, stdout = sys.stdout)
  
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


def extract_trees(datadir, results_family_dir, run_name, subst_model):
  results_dir = os.path.join(results_family_dir, "results")
  for family in fam.get_families_list(datadir):
    source = os.path.join(results_dir, family, "geneTree.newick")
    dest = fam.build_gene_tree_path_from_run(datadir, family, run_name)
    try:
      shutil.copy(source, dest)
    except:
      pass

def run(dataset, subst_model, strategy, starting_tree, cores, additional_arguments, resultsdir, do_analyze = True, do_extract = True):
  run_name = exp.getAndDelete("--run", additional_arguments, "generax-last." +subst_model) 
  arg_analyze = exp.getAndDelete("--analyze", additional_arguments, "yes")
  do_analyze = do_analyze and (arg_analyze == "yes") and (strategy != "EVAL")
  print("Run name " + run_name)
  sys.stdout.flush()
  mode = get_mode_from_additional_arguments(additional_arguments)
  if (not dataset in datasets):
    print("Error: " + dataset + " is not in " + str(datasets))
    exit(1)
  datadir = datasets[dataset]
  generax_families_file = os.path.join(resultsdir, "families.txt")
  build_generax_families_file(datadir, starting_tree, subst_model, generax_families_file)
  start = time.time()
  run_generax(datadir, strategy, generax_families_file, mode, cores, additional_arguments, resultsdir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "seqtimes") 
  if (do_extract):
    extract_trees(datadir, os.path.join(resultsdir, "generax"), run_name, subst_model)
  try:
    if (do_analyze):
      fast_rf_cells.analyze(datadir, "all", cores, run_name)
  except:
    print("Analyze failed!!!!")

  print("Output in " + resultsdir)

def launch(dataset, subst_model, strategy, starting_tree, cluster, cores, additional_arguments):
  command = ["python"]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("GeneRax", dataset, strategy + "_start_" + starting_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  submit_path = os.path.join(resultsdir, "submit.sh")
  command.append(resultsdir)
  exp.submit(submit_path, " ".join(command), cores, cluster) 
  

if (__name__ == "__main__"): 
  print("launch_generax " + str(sys.argv))
  is_run = ("--exprun" in sys.argv)
  resultsdir = ""
  if (is_run):
    resultsdir = sys.argv[-1]
    sys.argv = sys.argv[:-2]
    
  min_args_number = 7
  if (len(sys.argv) < min_args_number):
    for dataset in datasets:
      print("\t" + dataset)
    print("strategy: " + ",".join(get_possible_strategies()))
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset subst_model strategy starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
    sys.exit(1)

  dataset = sys.argv[1]
  subst_model = sys.argv[2]
  strategy = sys.argv[3]
  starting_tree = sys.argv[4]
  cluster = sys.argv[5]
  cores = int(sys.argv[6])
  additional_arguments = sys.argv[min_args_number:]
  check_inputs(strategy)

  if (starting_tree == "raxml"):
    print("use raxml-ng instead of raxml please")
    exit(1)

  if (is_run):
    run(dataset, subst_model, strategy, starting_tree, cores, additional_arguments, resultsdir)
  else:
    launch(dataset, subst_model, strategy, starting_tree, cluster, cores, additional_arguments)




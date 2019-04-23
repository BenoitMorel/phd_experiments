import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import saved_metrics
import experiments as exp
import shutil
import time
import fam

def get_speciesrax_datasets():
  root_datadir = os.path.join(exp.benoit_datasets_root, "families")
  datasets = {}
  for dataset in os.listdir(root_datadir):
    datasets[dataset] = os.path.join(root_datadir, dataset)
  return datasets

datasets = get_speciesrax_datasets()

def build_speciesrax_families_file(dataset, starting_tree, is_protein, output):
  families_dir = os.path.join(dataset, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    plop = 0
    for family in os.listdir(families_dir):
      
      family_path = os.path.join(families_dir, family)
      writer.write("- " + family + "\n")
      writer.write("starting_gene_tree = " + fam.get_gene_tree(family_path, starting_tree) + "\n")
      writer.write("alignment = " + fam.get_alignment_file(family_path) + "\n")
      writer.write("mapping = " + fam.get_mappings(dataset) + "\n")
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

def get_speciesrax_command(speciesrax_families_file, species_tree, additional_arguments, output_dir, mode, cores):
    executable = exp.speciesrax_exec
    if (mode == "gprof"):
      executable = exp.speciesrax_gprof_exec
    elif (mode == "scalasca"):
      executable = exp.speciesrax_scalasca_exec
    speciesrax_output = os.path.join(output_dir, "speciesrax")
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-f")
    command.append(speciesrax_families_file)
    command.append("-s")
    command.append(species_tree)
    command.append("-p")
    command.append(speciesrax_output)
    command.extend(additional_arguments)
    return " ".join(command)

def run_speciesrax(datadir, speciesrax_families_file, mode, cores, additional_arguments, resultsdir):
  species_tree = os.path.join(datadir, "speciesTree.newick")
  command = get_speciesrax_command(speciesrax_families_file, species_tree, additional_arguments, resultsdir, mode, cores)
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
    dest = os.path.join(output_msa_dir, prefix + ".newick")
    shutil.copy(source, dest)


def run(dataset, starting_tree, cores, additional_arguments, resultsdir):
  is_protein = exp.checkAndDelete("--protein", additional_arguments)
  run_name = exp.getAndDelete("--run", additional_arguments, "lastRun") 
  mode = get_mode_from_additional_arguments(additional_arguments)
  if (not dataset in datasets):
    print("Error: " + dataset + " is not in " + str(datasets))
    exit(1)
  datadir = datasets[dataset]
  speciesrax_families_file = os.path.join(resultsdir, "speciesrax_families.txt")
  build_speciesrax_families_file(datadir, starting_tree, is_protein, speciesrax_families_file)
  start = time.time()
  run_speciesrax(datadir, speciesrax_families_file, mode, cores, additional_arguments, resultsdir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 

  print("Output in " + resultsdir)

def launch(dataset, starting_tree, cluster, cores, additional_arguments):
  command = ["python"]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("SpeciesRax", dataset, "start_" + starting_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  submit_path = os.path.join(resultsdir, "submit.sh")
  command.append(resultsdir)
  exp.submit(submit_path, " ".join(command), cores, cluster) 
  

if (__name__ == "__main__"): 
  is_run = ("--exprun" in sys.argv)
  resultsdir = ""
  if (is_run):
    resultsdir = sys.argv[-1]
    sys.argv = sys.argv[:-2]
    
  min_args_number = 5
  if (len(sys.argv) < min_args_number):
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
    for dataset in datasets:
      print("\t" + dataset)
    print("starting_tree: " + ",".join(fam.get_possible_gene_trees()))
    sys.exit(1)

  dataset = sys.argv[1]
  starting_tree = sys.argv[2]
  cluster = sys.argv[3]
  cores = int(sys.argv[4])
  additional_arguments = sys.argv[min_args_number:]

  if (starting_tree == "raxml"):
    print("use raxml-ng instead of raxml please")
    exit(1)

  if (is_run):
    run(dataset, starting_tree, cores, additional_arguments, resultsdir)
  else:
    launch(dataset, starting_tree, cluster, cores, additional_arguments)




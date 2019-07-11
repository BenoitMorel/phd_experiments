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
from ete3 import Tree

def get_speciesrax_datasets():
  root_datadir = os.path.join(exp.benoit_datasets_root, "families")
  datasets = {}
  for dataset in os.listdir(root_datadir):
    datasets[dataset] = os.path.join(root_datadir, dataset)
  return datasets

datasets = get_speciesrax_datasets()

def build_speciesrax_families_file(dataset, starting_tree, subst_model, output):
  families_dir = os.path.join(dataset, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    plop = 0
    for family in os.listdir(families_dir):
      
      family_path = os.path.join(families_dir, family)
      writer.write("- " + family + "\n")
      writer.write("starting_gene_tree = " + fam.get_gene_tree(dataset, subst_model, family, starting_tree) + "\n")
      writer.write("alignment = " + fam.get_alignment_file(family_path) + "\n")
      writer.write("mapping = " + fam.get_mappings(dataset, family) + "\n")
      raxml_model = ""
      if (starting_tree != "random"):
        raxml_model = fam.get_raxml_best_model(dataset, subst_model, family)
      if (os.path.isfile(raxml_model)):
        writer.write("subst_model = " + raxml_model + "\n")
      else:
        writer.write("subst_model = " + subst_model + "\n")

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
  species_tree = fam.get_species_tree(datadir)
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

def av_rf(rf_cell):
  return float(rf_cell[0]) / float(rf_cell[1])

def analyze_results(datadir, resultsdir):
  true_species_tree = Tree(fam.get_species_tree(datadir), format = 1)
  starting_species_tree = Tree(os.path.join(resultsdir, "speciesrax", "starting_species_tree.newick"), format = 1)
  inferred_species_tree = Tree(os.path.join(resultsdir, "speciesrax", "inferred_species_tree.newick"), format = 1)
  starting_rooted_rf = true_species_tree.robinson_foulds(starting_species_tree, unrooted_trees = False)
  starting_unrooted_rf = true_species_tree.robinson_foulds(starting_species_tree, unrooted_trees = True)
  inferred_rooted_rf = true_species_tree.robinson_foulds(inferred_species_tree, unrooted_trees = False)
  inferred_unrooted_rf = true_species_tree.robinson_foulds(inferred_species_tree, unrooted_trees = True)
  print("Starting species tree:")
  print("  Rooted RF: " + str(av_rf(starting_rooted_rf)))
  print("  Unrooted RF: " + str(av_rf(starting_unrooted_rf)))
  print("Inferred species tree:")
  print("  Rooted RF: " + str(av_rf(inferred_rooted_rf)))
  print("  Unrooted RF: " + str(av_rf(inferred_unrooted_rf)))

def extract_results(datadir, subst_model, resultsdir, run_name):
  src = os.path.join(resultsdir, "speciesrax", "inferred_species_tree.newick")
  dest = fam.get_species_tree(datadir, subst_model, run_name)
  shutil.copyfile(src, dest)

def run(dataset, subst_model, starting_tree, cores, additional_arguments, resultsdir):
  print("Output in " + resultsdir)
  run_name = exp.getAndDelete("--run", additional_arguments, "speciesRax") 
  mode = get_mode_from_additional_arguments(additional_arguments)
  if (not dataset in datasets):
    print("Error: " + dataset + " is not in " + str(datasets))
    exit(1)
  datadir = datasets[dataset]
  speciesrax_families_file = os.path.join(resultsdir, "speciesrax_families.txt")
  build_speciesrax_families_file(datadir, starting_tree, subst_model, speciesrax_families_file)
  start = time.time()
  run_speciesrax(datadir, speciesrax_families_file, mode, cores, additional_arguments, resultsdir)
  saved_metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
  analyze_results(datadir, resultsdir) 
  extract_results(datadir, subst_model, resultsdir, run_name)

def launch(dataset, subst_model, starting_tree, cluster, cores, additional_arguments):
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
    
  min_args_number = 6
  if (len(sys.argv) < min_args_number):
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
    for dataset in datasets:
      print("\t" + dataset)
    sys.exit(1)

  dataset = sys.argv[1]
  subst_model = sys.argv[2]
  starting_tree = sys.argv[3]
  cluster = sys.argv[4]
  cores = int(sys.argv[5])
  additional_arguments = sys.argv[min_args_number:]

  if (starting_tree == "raxml"):
    print("use raxml-ng instead of raxml please")
    exit(1)

  if (is_run):
    run(dataset, subst_model, starting_tree, cores, additional_arguments, resultsdir)
  else:
    launch(dataset, subst_model, starting_tree, cluster, cores, additional_arguments)




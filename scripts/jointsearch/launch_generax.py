import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import analyze_dataset
import experiments as exp
import exp_jointsearch_utils as utils
import shutil

datasets = utils.get_generax_datasets()

def build_generax_families_file(dataset, starting_tree, is_protein, output):
  families_dir = os.path.join(dataset, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    plop = 0
    for family in os.listdir(families_dir):
      
      family_path = os.path.join(families_dir, family)
      writer.write("- " + family + "\n")
      writer.write("G: " + utils.get_gene_tree(family_path, starting_tree) + "\n")
      writer.write("A: " + utils.get_alignment_file(family_path) + "\n")
      writer.write("M: " + utils.get_mapping_file(family_path) + "\n")
      raxml_model = ""
      if (starting_tree != "random"):
        raxml_model = utils.get_raxml_model(family_path)
      if (os.path.isfile(raxml_model)):
        writer.write("L: " + raxml_model + "\n")
      else:
        if (is_protein):
          writer.write("L: LG\n")
        else:
          writer.write("L: GTR\n")

def run_generax(datadir, strategy, generax_families_file, mode, cores, additional_arguments, resultsdir):
  species_tree = os.path.join(datadir, "speciesTree.newick")
  command = utils.get_generax_command(generax_families_file, species_tree, strategy, additional_arguments, resultsdir, mode, cores)
  print(command)
  subprocess.check_call(command.split(" "))


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
  run_generax(datadir, strategy, generax_families_file, mode, cores, additional_arguments, resultsdir)
  extract_trees(os.path.join(datadir, "families"), os.path.join(resultsdir, "generax"), run_name)
  analyze_dataset.analyze(os.path.join(datadir, "families"), run_name)
  print("Output in " + resultsdir)

def launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments):
  command = ["python"]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("GeneRax", dataset, strategy + "_start_" + starting_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  result_msg = "GeneRax git: \n" + exp.get_git_info(exp.joint_search_root)
  exp.write_results_info(resultsdir, result_msg) 
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
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset strategy starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
    for dataset in datasets:
      print("\t" + dataset)
    print("strategy: " + ",".join(utils.get_possible_strategies()))
    print("starting_tree: " + ",".join(utils.get_possible_gene_trees()))
    sys.exit(1)

  dataset = sys.argv[1]
  strategy = sys.argv[2]
  starting_tree = sys.argv[3]
  cluster = sys.argv[4]
  cores = int(sys.argv[5])
  additional_arguments = sys.argv[min_args_number:]
  utils.check_inputs(starting_tree, strategy)

  if (is_run):
    run(dataset, strategy, starting_tree, cores, additional_arguments, resultsdir)
  else:
    launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments)




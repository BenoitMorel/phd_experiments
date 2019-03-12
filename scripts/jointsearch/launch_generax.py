import sys
import os
sys.path.insert(0, 'scripts')
import experiments as exp
import exp_jointsearch_utils as utils

datasets = utils.get_generax_datasets()

def build_generax_families_file(dataset, starting_tree, output):
  families_dir = os.path.join(dataset, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    for family in os.listdir(families_dir):
      family_path = os.path.join(families_dir, family)
      writer.write(family + "\n")
      writer.write("G: " + utils.get_gene_tree(family_path, starting_tree) + "\n")
      writer.write("A: " + utils.get_alignment_file(family_path) + "\n")
      writer.write("M: " + utils.get_mapping_file(family_path) + "\n")

def run_generax(dataset, strategy, generax_families_file, cores, additional_arguments, resultsdir):
  pass 

def run(dataset, strategy, starting_tree, cores, additional_arguments, resultsdir):
  datadir = datasets[dataset]
  generax_families_file = os.path.join(resultsdir, "generax_families.txt")
  build_generax_families_file(datadir, starting_tree, generax_families_file)
  run_generax(dataset, strategy, generax_families_file, cores, additional_arguments, resultsdir)

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




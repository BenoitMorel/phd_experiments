import sys
import os
import subprocess
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/msa_edition')
import experiments as exp
import exp_jointsearch_utils as utils
import msa_analyzer

datasets = utils.get_jointsearch_datasets()

def run_and_get_time(command):
  start = time.time()
  DEVNULL = open(os.devnull, 'wb', 0)
  subprocess.check_call(command, stdout=DEVNULL)
  elapsed = time.time() - start
  return elapsed


def run(dataset, strategy, starting_tree, cluster, cores, additional_arguments, resultsdir):
  datadir = datasets[dataset]
  gene_tree = utils.get_gene_tree(datadir, starting_tree)
  species_tree = os.path.join(datadir, "speciesTree.newick")
  alignment = os.path.join(datadir, "alignment.msa")
  mapping = os.path.join(datadir, "mapping.link")
  tree_size = msa_analyzer.get_taxa_number(alignment)
  print("Taxa number: " + str(tree_size))
  runs_number = 5
  for i in range(1, runs_number + 1):
    c = (cores * i) // runs_number
    run_dir = os.path.join(resultsdir, "cores_" + str(c))
    os.makedirs(run_dir)
    command = utils.get_jointsearch_command(gene_tree, species_tree, mapping, alignment, strategy, c, run_dir, False, additional_arguments).split(" ")
    t = run_and_get_time(command)
    print(str(c) + " cores: " + str(t))

def launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments):
  command = ["python"]
  command.extend(sys.argv)
  command.append("--exprun")
  resultsdir = os.path.join("JointSearchScaling", dataset, strategy + "_start_" + starting_tree, "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  result_msg = "JointSearch git: \n" + exp.get_git_info(exp.joint_search_root)
  exp.write_results_info(resultsdir, result_msg) 
  submit_path = os.path.join(resultsdir, "submit.sh")
  command.append(resultsdir)
  exp.submit(submit_path, " ".join(command), cores, cluster) 

is_run = ("--exprun" in sys.argv)
resultsdir = ""
if (is_run):
  resultsdir = sys.argv[-1]
  sys.argv = sys.argv[:-2]
  
min_args_number = 6
if (len(sys.argv) < min_args_number):
  print("Syntax error: python plot_jointsearch_scaling.py dataset strategy starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
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
  run(dataset, strategy, starting_tree, cluster, cores, additional_arguments, resultsdir)
else:
  launch(dataset, strategy, starting_tree, cluster, cores, additional_arguments)



import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp
import exp_jointsearch_utils as utils



def get_compare_rf_command(data_dir, tree, name, output):
  command = "echo Comparing true tree and " + name + " tree: >> " + output + "\n"
  command += "python " + exp.rf_distance_tool + " "
  command += utils.get_gene_tree(data_dir, "true") + " "
  command += utils.get_gene_tree(data_dir, tree) + " >> " + output
  return command



def launch_jointsearch(is_gprof):
  datasets = utils.get_jointsearch_datasets()
  max_args_number = 6
  if len(sys.argv) < max_args_number:
    print("Syntax error: python lauch_joint_search.py dataset strategy starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
    for dataset in datasets:
      print("\t" + dataset)
    print("strategy: " + ",".join(utils.get_possible_strategies()))
    print("starting_tree: " + ",".join(utils.get_possible_gene_trees()))
    sys.exit(0)

  dataset = sys.argv[1]
  strategy = sys.argv[2]
  starting_tree = sys.argv[3]
  cluster = sys.argv[4]
  cores = int(sys.argv[5])
  additional_arguments = sys.argv[max_args_number:]

  if (not (strategy in ["SPR", "NNI", "HYBRID"])):
    print("Unknown search strategy " + strategy)

  resultsdir = os.path.join("JointSearch", dataset, strategy + "_start_" + starting_tree, cluster + "_" + str(cores), "run")
  resultsdir = exp.create_result_dir(resultsdir, additional_arguments)
  result_msg = "JointSearch git: \n" + exp.get_git_info(exp.joint_search_root)
  exp.write_results_info(resultsdir, result_msg) 

  datadir = datasets[dataset]
  gene_tree = utils.get_gene_tree(datadir, starting_tree)
  species_tree = os.path.join(datadir, "speciesTree.newick")
  alignment = os.path.join(datadir, "alignment.msa")
  mapping = os.path.join(datadir, "mapping.link")
  output_dir = resultsdir 

  command = utils.get_jointsearch_command(gene_tree, species_tree, mapping, alignment, strategy, cores, output_dir, is_gprof, additional_arguments)

  result_tree = os.path.join(output_dir, "join_search.newick") 
  summary = os.path.join(output_dir, "summary.txt")


  command += "\n"
  command += "echo\n"
  command += get_compare_rf_command(datadir, "raxml", "raxml", summary)
  command += "\n"
  command += get_compare_rf_command(datadir, "treerecs", "treerecs", summary)
  command += "\n"
  command += get_compare_rf_command(datadir, result_tree, strategy + "_search", summary)
  command += "\n"
  command += "cat " + summary


  submit_path = os.path.join(resultsdir, "jointsearch_submit.sh")
  print("Results will be in " + resultsdir)
  exp.submit(submit_path, command, cores, cluster) 


if (__name__ == "__main__"): 
  launch_jointsearch(False)



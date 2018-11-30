import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp


def get_gene_tree(datadir, tree):
  if (tree == "raxml"):
    return os.path.join(datadir, "raxmlGeneTree.newick")
  elif (tree == "raxmls"):
    return os.path.join(datadir, "raxmlGeneTrees.newick")
  elif (tree == "true"):
    return os.path.join(datadir, "trueGeneTree.newick")
  elif (tree == "treerecs"):
    return os.path.join(datadir, "treerecsGeneTree.newick")
  elif (tree == "random"):
    return os.path.join(datadir, "randomGeneTree.newick")
  else:
    return tree


def get_tree_search_command(gene_tree, species_tree, mapping, alignment, strategy, cores, output_dir):
    executable = exp.joint_search_exec
    joint_search_output = os.path.join(output_dir, "join_search")
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-g")
    command.append(gene_tree)
    command.append("-a")
    command.append(alignment)
    command.append("-s")
    command.append(species_tree)
    command.append("-m")
    command.append(mapping)
    command.append("--strategy")
    command.append(strategy)
    command.append("-p")
    command.append(joint_search_output)
    return " ".join(command)

def get_compare_rf_command(data_dir, tree, name, output):
  command = "echo Comparing true tree and " + name + " tree: >> " + output + "\n"
  command += "python " + exp.rf_distance_tool + " "
  command += get_gene_tree(data_dir, "true") + " "
  command += get_gene_tree(data_dir, tree) + " >> " + output
  #command += "\nsleep 1"
  return command

root_datadir = os.path.join(exp.datasets_root, "joint_search")

datasets = {}
for dataset in os.listdir(root_datadir):
  datasets[dataset] = os.path.join(root_datadir, dataset)

max_args_number = 6
if len(sys.argv) < max_args_number:
  print("Syntax error: python lauch_joint_search.py dataset strategy starting_tree cluster cores [additional paremeters].\n Suggestions of datasets: ")
  for dataset in datasets:
    print("\t" + dataset)
  print("Cluster can be either normal, haswell or magny")
  print("strategy: SPR, NNI, HYBRID")
  print("starting_tree: raxml, raxmls, true, treerecs, random")
  sys.exit(0)

dataset = sys.argv[1]
strategy = sys.argv[2]
starting_tree = sys.argv[3]
cluster = sys.argv[4]
cores = int(sys.argv[5])

if (not (strategy in ["SPR", "NNI", "HYBRID"])):
  print("Unknown search strategy " + strategy)

resultsdir = os.path.join("JointSearch", dataset, strategy + "_start_" + starting_tree, cluster + "_" + str(cores), "run")
for i in range(max_args_number, len(sys.argv)):
  resultsdir += "_" + sys.argv[i]

resultsdir = exp.create_result_dir(resultsdir)
result_msg = "JointSearch git: \n" + exp.get_git_info(exp.joint_search_root)
exp.write_results_info(resultsdir, result_msg) 

datadir = datasets[dataset]
gene_tree = get_gene_tree(datadir, starting_tree)
species_tree = os.path.join(datadir, "speciesTree.newick")
alignment = os.path.join(datadir, "alignment.msa")
mapping = os.path.join(datadir, "phyldog", "phyldogMapping.link")
#smap = os.path.join(datadir, "mapping.txt")
#if (not os.path.isfile(smap)):
#  smap = ""
output_dir = resultsdir 

command = get_tree_search_command(gene_tree, species_tree, mapping, alignment, strategy, cores, output_dir)
for i in range(max_args_number, len(sys.argv)):
  command += " " + sys.argv[i]

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




import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp


def get_gene_tree(datadir, tree):
  if (tree == "raxml"):
    return os.path.join(datadir, "raxmlGeneTree.newick")
  elif (tree == "true"):
    return os.path.join(datadir, "trueGeneTree.newick")
  elif (tree == "treerecs"):
    return os.path.join(datadir, "treerecsGeneTree.newick")
  elif (tree == "random"):
    return os.path.join(datadir, "randomGeneTree.newick")
  else:
    return tree


def get_tree_search_command(gene_tree, species_tree, alignment, strategy, cores, output_dir):
    executable = exp.joint_search_exec
    joint_search_output = os.path.join(output_dir, "join_search")
    command = []
    command.append(executable)
    command.append(gene_tree)
    command.append(alignment)
    command.append(species_tree)
    command.append(strategy)
    command.append(str(cores))
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


if len(sys.argv) != 6:
  print("Syntax error: python lauch_treerecs.py dataset strategy starting_tree cluster cores.\n Suggestions of datasets: ")
  for dataset in datasets:
    print("\t" + dataset)
  print("Cluster can be either normal, haswell or magny")
  print("strategy: SPR or NNI")
  print("starting_tree: raxml, true, treerecs, random")
  sys.exit(0)

dataset = sys.argv[1]
strategy = sys.argv[2]
starting_tree = sys.argv[3]
cluster = sys.argv[4]
cores = int(sys.argv[5])

if (not (strategy in ["SPR", "NNI"])):
  print("Unknown search strategy " + strategy)

resultsdir = exp.create_result_dir(os.path.join("treerecs", "joint_search", dataset, strategy + "_start_" + starting_tree, cluster + "_" + str(cores), "run"))
result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 

datadir = datasets[dataset]
gene_tree = get_gene_tree(datadir, starting_tree)
species_tree = os.path.join(datadir, "speciesTree.newick")
alignment = os.path.join(datadir, "alignment.msa")
#smap = os.path.join(datadir, "mapping.txt")
#if (not os.path.isfile(smap)):
#  smap = ""
output_dir = resultsdir 

command = get_tree_search_command(gene_tree, species_tree, alignment, strategy, cores, output_dir)

result_tree = command.split(" ")[-1] + ".newick"
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


submit_path = os.path.join(resultsdir, "treerecs_submit.sh")
print("Results will be in " + resultsdir)
exp.submit(submit_path, command, cores, cluster) 




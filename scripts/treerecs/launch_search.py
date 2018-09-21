import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp



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

root_datadir = os.path.join(exp.datasets_root, "joint_search")

datasets = {}
for dataset in os.listdir(root_datadir):
  datasets[dataset] = os.path.join(root_datadir, dataset)


if len(sys.argv) != 7:
  print("Syntax error: python lauch_treerecs.py dataset run_name strategy starting_tree cluster cores.\n Suggestions of datasets: ")
  for dataset in datasets:
    print("\t" + dataset)
  print("Cluster can be either normal, haswell or magny")
  print("strategy: SPR or NNI")
  print("starting_tree: raxml, true, treerecs, random")
  sys.exit(0)

dataset = sys.argv[1]
run_name = sys.argv[2]
strategy = sys.argv[3]
starting_tree = sys.argv[4]
cluster = sys.argv[5]
cores = int(sys.argv[6])

if (not (strategy in ["SPR", "NNI"])):
  print("Unknown search strategy " + strategy)

resultsdir = exp.create_result_dir(os.path.join("treerecs", "joint_search", dataset, run_name + "_" + strategy + "_start_" + starting_tree, cluster + "_" + str(cores), "run"))
result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 

datadir = datasets[dataset]
gene_tree = ""
if (starting_tree == "raxml"):
  gene_tree = os.path.join(datadir, "raxmlGeneTree.newick")
elif (starting_tree == "true"):
  gene_tree = os.path.join(datadir, "trueGeneTree.newick")
elif (starting_tree == "treerecs"):
  gene_tree = os.path.join(datadir, "treerecsGeneTree.newick")
elif (starting_tree == "random"):
  gene_tree = os.path.join(datadir, "randomGeneTree.newick")
else:
  print("Error: wrong starting tree")
  sys.exit(1)
species_tree = os.path.join(datadir, "speciesTree.newick")
alignment = os.path.join(datadir, "alignment.msa")
#smap = os.path.join(datadir, "mapping.txt")
#if (not os.path.isfile(smap)):
#  smap = ""
output_dir = resultsdir 

command = get_tree_search_command(gene_tree, species_tree, alignment, strategy, cores, output_dir)
submit_path = os.path.join(resultsdir, "treerecs_submit.sh")
print("Results will be in " + resultsdir)
exp.submit(submit_path, command, cores, cluster) 





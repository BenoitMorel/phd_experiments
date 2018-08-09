import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp



def get_treerecs_command(gene_tree, species_tree, alignment, smap, thresholds_number, tree_search, output_dir, cores):
    treerecs_exec = exp.treerecs_exec
    treerecs_output = os.path.join(output_dir, "treerecs_output")
    command = []
    command.append(treerecs_exec)
    command.append("-g")
    command.append(gene_tree)
    command.append("-s")
    command.append(species_tree)
    command.append("-o")
    command.append(treerecs_output)
    command.append("-a")
    command.append(alignment)
    command.append("-t")
    command.append("all")
    command.append("--ale-evaluation")
    command.append("-T")
    command.append(str(thresholds_number))
    command.append("--select-best-tree")
    if (tree_search):
      command.append("--tree-search")
    if (smap):
      command.append("-S")
      command.append(smap)
    if (cores != 1):
      command.append("-P")
      command.append(str(cores))
    return " ".join(command)

root_datadir = os.path.join(exp.datasets_root, "treerecs")

datasets = {}
for dataset in os.listdir(root_datadir):
  datasets[dataset] = os.path.join(root_datadir, dataset)
datasets["simuls"] = os.path.join(exp.bigdatasets_root, "simuls") 
datasets["simuls_higher_rate"] = os.path.join(exp.bigdatasets_root, "simuls_higher_rate") 
datasets["sub_simuls_higher_rate"] = os.path.join(exp.bigdatasets_root, "sub_simuls_higher_rate") 


if len(sys.argv) != 6:
  print("Syntax error: python lauch_treerecs.py dataset tree-search run_time cluster cores.\n Suggestions of datasets: ")
  for dataset in datasets:
    print("\t" + dataset)
  print("Cluster can be either normal, haswell or magny")
  sys.exit(0)

dataset = sys.argv[1]
tree_search = int(sys.argv[2])
run_name = sys.argv[3]
cluster = sys.argv[4]
cores = int(sys.argv[5])

subdir_name = dataset
if (tree_search != 0):
  subdir_name += "_treesearch"

resultsdir = exp.create_result_dir(os.path.join("treerecs", "launch_treerecs",run_name, subdir_name, cluster + "_" + str(cores), "run"))
result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 

datadir = datasets[dataset]
gene_tree = os.path.join(datadir, "geneTrees.newick")
species_tree = os.path.join(datadir, "speciesTree.newick")
alignment = os.path.join(datadir, "alignment.txt")
smap = os.path.join(datadir, "mapping.txt")
if (not os.path.isfile(smap)):
  smap = ""
thresholds_number = 7
output_dir = resultsdir 





command = get_treerecs_command(gene_tree, species_tree, alignment, smap, thresholds_number, tree_search, output_dir, cores)
submit_path = os.path.join(resultsdir, "treerecs_submit.sh")
print("Results will be in " + resultsdir)
exp.submit(submit_path, command, cores, cluster) 





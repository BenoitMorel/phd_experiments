import sys
import os
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp



def get_evaluator_command(gene_tree, species_tree, alignment, smap, output_dir, cores):
    evaluator_output = os.path.join(output_dir, "evaluator_output.txt")
    command = []
    command.append(exp.joint_likelihood_evaluator_exec)
    command.append("-g")
    command.append(gene_tree)
    command.append("-s")
    command.append(species_tree)
    command.append("-o")
    command.append(evaluator_output)
    command.append("-a")
    command.append(alignment)
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
  print("Syntax error: python lauch_joint_ll_evaluator.py dataset geneTree run_name cluster cores.\n Suggestions of datasets: ")
  for dataset in datasets:
    print("\t" + dataset)
  print("Cluster can be either normal, haswell or magny")
  sys.exit(0)

dataset = sys.argv[1]
gene_tree = sys.argv[2]
run_name = sys.argv[3]
cluster = sys.argv[4]
cores = int(sys.argv[5])


resultsdir = exp.create_result_dir(os.path.join("treerecs", "launch_joint_likelihood_evaluator", run_name, dataset, cluster + "_" + str(cores), "run"))
result_msg = "Treerecs git: \n" + exp.get_git_info(exp.treerecs_root)
exp.write_results_info(resultsdir, result_msg) 

datadir = datasets[dataset]
species_tree = os.path.join(datadir, "speciesTree.newick")
alignment = os.path.join(datadir, "alignment.txt")
smap = os.path.join(datadir, "mapping.txt")
if (not os.path.isfile(smap)):
  smap = ""





command = get_evaluator_command(gene_tree, species_tree, alignment, smap,  resultsdir, abs(cores))
submit_path = os.path.join(resultsdir, "joint_submit.sh")
print("Results will be in " + resultsdir)
exp.submit(submit_path, command, cores, cluster) 






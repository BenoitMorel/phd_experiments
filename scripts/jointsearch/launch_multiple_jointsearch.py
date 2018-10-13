import sys
import os
sys.path.insert(0, 'scripts')
import experiments as exp

#todobenoit this is ugly
def get_nodes_number(gene_tree_file):
  lines = open(gene_tree_file).readlines()
  line = lines[0]
  return line.count("(")


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

def generate_scheduler_commands_file(families_dir, starting_tree, strategy, nodes_per_core, additional_arguments, cores, output_dir, scheduler_output_dir):
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  results_dir = os.path.join(scheduler_output_dir, "results")
  with open(scheduler_commands_file, "w") as writer:
    for family in os.listdir(families_dir):
      family_path = os.path.join(families_dir, family)
      gene_tree = get_gene_tree(family_path, starting_tree)
      if (len(open(gene_tree).readlines()) == 0):
        continue
      tree_size = get_nodes_number(gene_tree)
      if (0 == nodes_per_core):
        writer.write(family + " 1 " + str(tree_size) + " ")
      else:
        writer.write(family + " " + str(max(1, int(tree_size / nodes_per_core))) + " " + str(tree_size) + " ")
      writer.write("-g " + gene_tree + " ")
      writer.write("-s " + os.path.join(family_path, "speciesTree.newick") +  " ")
      writer.write("-a " + os.path.join(family_path, "alignment.msa") +  " ")
      os.makedirs(os.path.join(results_dir, family))
      writer.write("-p " + os.path.join(results_dir, family, "jointsearch") + " ") 
      writer.write("--strategy " + strategy + " ")
      writer.write(additional_arguments + " ")
      writer.write("\n")

  return scheduler_commands_file


def generate_scheduler_command(command_file, parallelization, cores, scheduler_output_dir):
  command = ""
  isMPI = (parallelization == "onecore") or (parallelization == "split")
  if (isMPI):
    command += "mpirun -np " + str(cores) + " "
  command += exp.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  if (parallelization == "split"):
    command += exp.joint_search_lib + " "
  elif (parallelization == "openmp" or parallelization == "onecore"):
    command += exp.joint_search_exec + " "
  else:
    print("Unknown scheduler parallelization " + parallelization)
    sys.exit(0)
  command += command_file + " "
  command += scheduler_output_dir + " " 
  command += "0 "
  return command 

def generate_analyse_command(families_dir, pargenes_dir):
  command = "python "
  command += exp.analyse_pargenes_jointsearch_tool + " "
  command += families_dir + " " +  pargenes_dir
  return command

max_args_number = 6
if len(sys.argv) < max_args_number:
  print("Syntax error: python launch_multiple_jointsearch.py strategy starting_tree scheduler_implem cluster cores [additional paremeters].")
  print("Cluster can be either normal, haswell or magny")
  print("strategy: SPR, NNI, HYBRID")
  print("starting_tree: raxml, true, treerecs, random")
  print("scheduler implem: split, onecore, openmp")
  sys.exit(0)


strategy = sys.argv[1]
starting_tree = sys.argv[2]
parallelization = sys.argv[3]
cluster = sys.argv[4]
cores = int(sys.argv[5])
nodes_per_core = 20

if (not (strategy in ["SPR", "NNI", "HYBRID"])):
  print("Unknown search strategy " + strategy)

resultsdir = os.path.join("MultipleJointSearch", strategy + "_" + str(nodes_per_core) + "_start_" + starting_tree + "_" + parallelization, cluster + "_" + str(cores), "run")
for i in range(max_args_number, len(sys.argv)):
  resultsdir += "_" + sys.argv[i]

resultsdir = exp.create_result_dir(resultsdir)
result_msg = "JointSearch git: \n" + exp.get_git_info(exp.joint_search_root)
exp.write_results_info(resultsdir, result_msg) 
output_dir = resultsdir 

datadir = os.path.join(exp.bigdatasets_root, "simuls_higher_rate")
families_dir = os.path.join(datadir, "families")


additional_arguments = ""
for i in range(max_args_number, len(sys.argv)):
  additional_arguments += " " + sys.argv[i]
scheduler_output_dir = os.path.join(output_dir, "scheduler_run")
os.makedirs(scheduler_output_dir)
scheduler_commands_file = generate_scheduler_commands_file(families_dir, starting_tree, strategy, nodes_per_core, additional_arguments, cores, output_dir, scheduler_output_dir)
command = generate_scheduler_command(scheduler_commands_file, parallelization, cores,  scheduler_output_dir)

command += "\n" + generate_analyse_command(families_dir, scheduler_output_dir)

submit_path = os.path.join(resultsdir, "jointsearch_scheduler_submit.sh")
print("Results will be in " + resultsdir)
exp.submit(submit_path, command, cores, cluster) 





import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp


def get_possible_strategies():
  return ["SPR", "EVAL"]



def get_jointsearch_datasets():
  root_datadir = os.path.join(exp.datasets_root, "joint_search")
  datasets = {}
  for dataset in os.listdir(root_datadir):
    datasets[dataset] = os.path.join(root_datadir, dataset)
  return datasets


def get_jointsearch_command(gene_tree, species_tree, mapping, alignment, strategy, cores, output_dir, mode, additional_arguments):
    executable = exp.joint_search_exec
    if (mode == "gprof"):
      executable = exp.joint_search_gprof_exec
    elif (mode == "scalasca"):
      executable = exp.joint_search_scalasca_exec
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
    command.extend(additional_arguments)
    return " ".join(command)




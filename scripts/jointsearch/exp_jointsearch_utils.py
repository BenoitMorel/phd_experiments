import os
import sys
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
    return "__random__";
  else:
    return tree


def get_jointsearch_command(gene_tree, species_tree, mapping, alignment, strategy, cores, output_dir, is_gprof, additional_arguments):
    executable = exp.joint_search_exec
    if (is_gprof):
      executable = exp.joint_search_gprof_exec
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


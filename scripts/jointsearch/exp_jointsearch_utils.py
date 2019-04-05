import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp

def get_possible_gene_trees():
  return ["raxml", "raxmls", "true", "treerecs", "notung", "phyldog", "random", "ALE-D(T)L", "GeneRax-D(T)L-[Random, Raxml]"]

def get_possible_strategies():
  return ["SPR", "EVAL"]

def get_gene_tree(familydir, tree):
  gene_trees_dir = os.path.join(familydir, "gene_trees")
  lower_tree = tree.lower()
  if (lower_tree == "raxml-ng"):
    return os.path.join(gene_trees_dir, "raxmlGeneTree.newick")
  elif (lower_tree == "raxmls"):
    return os.path.join(gene_trees_dir, "raxmlGeneTrees.newick")
  elif (lower_tree == "true"):
    return os.path.join(familydir, "trueGeneTree.newick")
  elif (lower_tree == "treerecs"):
    return os.path.join(gene_trees_dir, "treerecsGeneTree.newick")
  elif (lower_tree == "phyldog"):
    return os.path.join(gene_trees_dir, "phyldogGeneTree.newick")
  elif (lower_tree == "notung"):
    return os.path.join(gene_trees_dir, "notungGeneTree.newick")
  elif ("GeneRax" in tree):
    return os.path.join(familydir, "results", tree + ".newick")
  elif ("ALE" in tree):
    return os.path.join(gene_trees_dir, tree + "GeneTree.newick")
  elif (tree == "random"):
    return "__random__";
  else:
    return tree

def get_mapping_file(datadir):
  return os.path.join(datadir, "mapping.link")

def get_alignment_file(datadir):
  return os.path.join(datadir, "alignment.msa")

def get_raxml_model(datadir):
  return os.path.join(datadir, "raxmlBestModel.txt")

def check_inputs(starting_tree, strategy):
  if (not (strategy in get_possible_strategies())):
    print("Unknown search strategy " + strategy)
    exit(1)

def get_jointsearch_datasets():
  root_datadir = os.path.join(exp.datasets_root, "joint_search")
  datasets = {}
  for dataset in os.listdir(root_datadir):
    datasets[dataset] = os.path.join(root_datadir, dataset)
  return datasets

def get_generax_datasets():
  root_datadir = os.path.join(exp.benoit_datasets_root, "families")
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

def get_generax_command(generax_families_file, species_tree, strategy, additional_arguments, output_dir, mode, cores):
    executable = exp.generax_exec
    if (mode == "gprof"):
      executable = exp.generax_gprof_exec
    elif (mode == "scalasca"):
      executable = exp.generax_scalasca_exec
    generax_output = os.path.join(output_dir, "generax")
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-f")
    command.append(generax_families_file)
    command.append("-s")
    command.append(species_tree)
    command.append("--strategy")
    command.append(strategy)
    command.append("-p")
    command.append(generax_output)
    command.extend(additional_arguments)
    return " ".join(command)



import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import perturbate

def mkdir(directory):
  try:
    os.makedirs(directory)
  except:
    pass

######################
# global directories
######################

def get_species_dir(datadir):
  return os.path.join(datadir, "species_trees")

def get_adjacencies_dir(datadir):
  return os.path.join(datadir, "adjacencies")

def get_families_dir(datadir):
  return os.path.join(datadir, "families")

def get_alignments_dir(datadir):
  return os.path.join(datadir, "alignments")

######################
# global files
#####################

def get_species_tree(datadir):
  return os.path.join(get_species_dir(datadir), "speciesTree.newick")

def get_true_species_tree(datadir):
  return os.path.join(get_species_dir(datadir), "trueSpeciesTree.newick")

def get_phyldog_species_tree(datadir):
  return os.path.join(get_species_dir(datadir), "phyldogSpeciesTree.newick")

def get_adjacencies(datadir):
  return os.path.join(get_adjacencies_dir(datadir), "adjacencies.txt")

def get_prefixed_adjacencies(datadir):
  return get_adjacencies(datadir) + ".prefixed"

def get_deco_mappings(datadir):
  return os.path.join(get_adjacencies_dir(datadir), "deco_mappings.txt")


#####################
# families
#####################

def get_families_list(datadir):
  return os.listdir(get_families_dir(datadir))

def get_family_path(datadir, family):
  return os.path.join(get_families_dir(datadir), family)

#####################
# per-family directories
#####################

def get_gene_tree_dir(datadir, family):
  return os.path.join(get_family_path(datadir, family), "gene_trees")

def get_misc_dir(datadir, family):
  return os.path.join(get_family_path(datadir, family), "misc")

def get_mappings_dir(datadir, family):
  return os.path.join(get_family_path(datadir, family), "mappings")

#####################
# per-family files
#####################

def get_alignment(datadir, family):
  return os.path.join(get_family_path(datadir, family), "alignment.msa")

def get_true_tree(datadir, family):
  return os.path.join(get_gene_tree_dir(datadir, family), "trueGeneTree.newick")

def get_random_tree(datadir, family):
  return os.path.join(get_gene_tree_dir(datadir, family), "randomGeneTree.newick")

def get_raxml_tree(datadir, family):
  return os.path.join(get_gene_tree_dir(datadir, family), "raxmlGeneTree.newick")

def get_raxml_multiple_trees(datadir, family):
  return os.path.join(get_gene_tree_dir(datadir, family), "raxmlGeneTrees.newick")

def get_phyldog_tree(datadir, family):
  return os.path.join(get_gene_tree_dir(datadir, family), "phyldogGeneTree.newick")

def get_treerecs_tree(datadir, family):
  return os.path.join(get_gene_tree_dir(datadir, family), "treerecsGeneTree.newick")

def get_notung_tree(datadir, family, threshold):
  return os.path.join(get_gene_tree_dir(datadir, family), "notung" + str(threshold) + "GeneTree.newick")

def get_ale_tree(datadir, family, method):
  return os.path.join(get_gene_tree_dir(datadir, family), method + "GeneTree.newick")

def get_mappings(datadir, family):
  return os.path.join(get_mappings_dir(datadir, family), "mapping.link")

def get_treerecs_mappings(datadir, family):
  return os.path.join(get_mappings_dir(datadir, family), "treerecs_mapping.link")

def get_alignment_file(datadir):
  return os.path.join(datadir, "alignment.msa")

def get_raxml_model(datadir):
  return os.path.join(datadir, "raxmlBestModel.txt")


def get_gene_tree(familydir, tree):
  gene_trees_dir = os.path.join(familydir, "gene_trees")
  lower_tree = tree.lower()
  if (lower_tree == "raxml-ng"):
    return os.path.join(gene_trees_dir, "raxmlGeneTree.newick")
  elif (lower_tree == "raxmls"):
    return os.path.join(gene_trees_dir, "raxmlGeneTrees.newick")
  elif (lower_tree == "true"):
    return os.path.join(gene_trees_dir, "trueGeneTree.newick")
  elif (lower_tree == "treerecs"):
    return os.path.join(gene_trees_dir, "treerecsGeneTree.newick")
  elif (lower_tree == "phyldog"):
    return os.path.join(gene_trees_dir, "phyldogGeneTree.newick")
  elif ("notung" in lower_tree):
    return os.path.join(gene_trees_dir, lower_tree + "GeneTree.newick")
  elif ("generax" in lower_tree):
    return os.path.join(familydir, "results", lower_tree + "GeneTree.newick")
  elif ("ale" in lower_tree):
    return os.path.join(gene_trees_dir, tree + "GeneTree.newick")
  elif (lower_tree == "random"):
    res = os.path.join(gene_trees_dir, "randomGeneTree.newick")
    if (os.path.isfile(res)):
        return res
    else:
      return "__random__";
  else:
    return tree

def get_method_name_aux(tree_path):
  tree = os.path.basename(tree_path)
  if (tree == "raxmlGeneTree.newick"):
    return "raxml-ng"
  elif (tree == "raxmlGeneTrees.newick"):
    return "raxmls"
  elif (tree == "trueGeneTree.newick"):
    return "true"
  elif (tree == "treerecsGeneTree.newick"):
    return "treerecs"
  elif (tree == "phyldogGeneTree.newick"):
    return "phyldog"
  elif (tree == "randomGeneTree.newick"):
    return "random"
  elif (tree == "randomGeneTree.newick"):
    return "random"
  else:
    return tree.replace("GeneTree.newick", "").replace(".newick", "")
  

def get_method_name(family_dir, tree_path):
  if (not tree_path.endswith(".newick")):
    return None
  method_name = get_method_name_aux(tree_path)
  if (not os.path.isfile(get_gene_tree(family_dir, method_name))):
    return None
  return method_name

def get_ran_methods(datadir):
  methods = []
  families_dir = os.path.join(datadir, "families")
  one_family_dir = os.path.join(families_dir, os.listdir(families_dir)[0])
  print(one_family_dir)
  directories_to_check = []
  directories_to_check.append(one_family_dir)
  directories_to_check.append(os.path.join(one_family_dir, "gene_trees"))
  directories_to_check.append(os.path.join(one_family_dir, "results"))

  for directory in directories_to_check:
    for tree in os.listdir(directory):
      method = get_method_name(one_family_dir, os.path.join(directory, tree))
      if (method != None and method != "raxmls"):
        methods.append(method)
  return methods

def get_param_from_dataset_name(parameter, dataset):
  if (parameter == "species"):
    return dataset.split("_")[1][1:]
  if (parameter == "families"):
    return dataset.split("_")[2][1:]
  elif (parameter == "sites"):
    return dataset.split("_")[3][5:]
  elif (parameter == "model"):
    return dataset.split("_")[4]
  elif (parameter == "bl"):
    return dataset.split("_")[5][2:]
  elif (parameter == "dup_rate"):
    return dataset.split("_")[6][1:]
  elif (parameter == "loss_rate"):
    return dataset.split("_")[7][1:]
  elif (parameter == "transfer_rate"):
    return dataset.split("_")[8][1:]
  elif (parameter == "perturbation"):
    return dataset.split("_")[9][1:]
  elif (parameter == "tl_ratio"):
    t = get_param_from_dataset_name("transfer_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    return  str(float(t)/float(l))
  elif (parameter == "dl_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    return str(float(d)/float(l))
  elif (parameter == "dt_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    t = get_param_from_dataset_name("transfer_rate", dataset)
    return str(float(d)/float(t))
  elif (parameter == "av_rate"):
    d = float(get_param_from_dataset_name("dup_rate", dataset))
    l = float(get_param_from_dataset_name("loss_rate", dataset))
    t = float(get_param_from_dataset_name("transfer_rate", dataset))
    return str((d + t + l) / 2.0)
  else:
    return "invalid"

######################
# convert functions
######################
def convert_to_phyldog_species_tree(speciesTree, phyldog_species_tree):
  command = "sed s/)[nHR][0123456789a-zA-Z]*/)/g " + speciesTree 
  with open(phyldog_species_tree, "w") as output:
    subprocess.check_call(command.split(" "), stdout=output)
  print(open(speciesTree).read())
  subprocess.check_call(command.split(" "))
  print("")
  print("")

def convert_phyldog_to_treerecs_mapping(phyldog_mappings, treerecs_mappings):
  lines = open(phyldog_mappings).readlines()
  with open(treerecs_mappings, "w") as writer:
    for line in lines:
      split = line.split(":")
      species = split[0]
      genes = split[1].split(";")
      for gene in genes:
        writer.write(gene.replace("\n", "") + " " + species + "\n")

def write_phyldog_mapping(species_to_genes_dict, output_file):
  with open(output_file, "w") as writer:
    for species in species_to_genes_dict:
      writer.write(species + ":" + ";".join(species_to_genes_dict[species]) + "\n")

#######################
#  Directory helpers
######################

def init_top_directories(datadir):
  mkdir(datadir)
  mkdir(get_species_dir(datadir))
  mkdir(get_adjacencies_dir(datadir))
  mkdir(get_families_dir(datadir))
  mkdir(get_alignments_dir(datadir))

def init_family_directories(datadir, family):
  mkdir(get_gene_tree_dir(datadir, family))
  mkdir(get_misc_dir(datadir, family))
  mkdir(get_mappings_dir(datadir, family))

def init_families_directories(datadir, families):
  for family in families:
    init_family_directories(datadir, family)


def postprocess_datadir(datadir):
  # phyldog species trees
  convert_to_phyldog_species_tree(get_species_tree(datadir), get_phyldog_species_tree(datadir)) 
  # alignments
  for family in get_families_list(datadir):
    family_dir = get_family_path(datadir, family)
    exp.relative_symlink(get_alignment(datadir, family), os.path.join(datadir, "alignments", family + ".fasta"))
    convert_phyldog_to_treerecs_mapping(get_mappings(datadir, family), get_treerecs_mappings(datadir, family)) 

def perturbate_species_tree(datadir, perturbation):
  shutil.copyfile(get_species_tree(datadir), get_true_species_tree(datadir))
  perturbate.perturbate(get_true_species_tree(datadir), get_species_tree(datadir), perturbation)





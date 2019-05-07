import os
import sys
import subprocess

def mkdir(directory):
  try:
    os.makedirs(directory)
  except:
    pass

def getSpeciesTree(dataset_dir):
  return os.path.join(dataset_dir, "speciesTree.newick")

def get_adjacencies(dataset_dir):
  return os.path.join(dataset_dir, "adjacencies.txt")

def get_prefixed_adjacencies(dataset_dir):
  return get_adjacencies(dataset_dir) + ".prefixed"

def get_deco_mappings(dataset_dir):
  return os.path.join(dataset_dir, "deco_mappings.txt")

def getFamiliesDir(dataset_dir):
  return os.path.join(dataset_dir, "families")

def getFamiliesList(dataset_dir):
  return os.listdir(getFamiliesDir(dataset_dir))

def getFamily(dataset_dir, family):
  return os.path.join(getFamiliesDir(dataset_dir), family)

def getTreesDir(dataset_dir, family):
  return os.path.join(getFamily(dataset_dir, family), "gene_trees")

def getMiscDir(dataset_dir, family):
  return os.path.join(getFamily(dataset_dir, family), "misc")

def getAlignment(dataset_dir, family):
  return os.path.join(getFamily(dataset_dir, family), "alignment.msa")

def getMappings(dataset_dir, family):
  return os.path.join(getFamily(dataset_dir, family), "mappings.link")


def getTrueTree(dataset_dir, family):
  return os.path.join(getFamily(dataset_dir, family), "trueGeneTree.newick")

def getRandomTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "randomGeneTree.newick")

def getRaxmlTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "raxmlGeneTree.newick")

def getRaxmlMultipleTrees(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "raxmlGeneTrees.newick")

def getPhyldogTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "phyldogGeneTree.newick")

def getTreerecsTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "treerecsGeneTree.newick")

def getNotungTree(dataset_dir, family, threshold = None):
  if (threshold == None):
    return os.path.join(getTreesDir(dataset_dir, family), "notungGeneTree.newick")
  else:
    return os.path.join(getTreesDir(dataset_dir, family), "notung" + str(threshold) + "GeneTree.newick")

def getALETree(dataset_dir, family, method):
  return os.path.join(getTreesDir(dataset_dir, family), method + "GeneTree.newick")



def init_dataset_dir(dataset_dir):
  for family in getFamiliesList(dataset_dir):
    mkdir(getTreesDir(dataset_dir, family))
    mkdir(getMiscDir(dataset_dir, family))

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

def get_species_tree(datadir):
  return os.path.join(datadir, "speciesTree.newick")

def get_phyldog_species_tree(datadir):
  return os.path.join(datadir, "phyldogSpeciesTree.newick")

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
  elif ("notung" in lower_tree):
    return os.path.join(gene_trees_dir, lower_tree + "GeneTree.newick")
  elif ("generax" in lower_tree):
    return os.path.join(familydir, "results", tree + ".newick")
  elif ("rec" in tree):
    return os.path.join(familydir, "results", tree + ".newick")
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

def get_possible_gene_trees():
  return ["raxml", "raxmls", "true", "treerecs", "notung", "phyldog", "random", "ALE-D(T)L", "GeneRax-D(T)L-[Random, Raxml]"]

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

def get_successfully_ran_methods(datadir):
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
      if (method != None):
        methods.append(method)
  return methods

def get_possible_gene_trees():
  return ["raxml", "raxmls", "true", "treerecs", "notung", "phyldog", "random", "ALE-D(T)L", "GeneRax-D(T)L-[Random, Raxml]"]

def get_mappings(datadir, family):
  return os.path.join(getFamily(datadir, family), "mapping.link")

def get_treerecs_mappings(datadir, family):
  return os.path.join(getFamily(datadir, family), "treerecs_mapping.link")


def get_alignment_file(datadir):
  return os.path.join(datadir, "alignment.msa")

def get_raxml_model(datadir):
  return os.path.join(datadir, "raxmlBestModel.txt")

def convertToPhyldogSpeciesTree(speciesTree, phyldogSpeciesTree):
  command = "sed s/)[nHR][0123456789a-zA-Z]*/)/g " + speciesTree #+ " > " + phyldogSpeciesTree
  with open(phyldogSpeciesTree, "w") as output:
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


if (__name__ == "__main__"): 
  if (len(sys.argv) != 2): 
     print("Syntax: python " + os.path.basename(__file__) + " param1")
     exit(1)
  print(get_successfully_ran_methods(sys.argv[1]))





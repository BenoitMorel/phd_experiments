import os
import sys

def mkdir(directory):
  try:
    os.makedirs(directory)
  except:
    pass

def getTreesDir(dataset_dir, family):
  return os.path.join(dataset_dir, "families", family, "gene_trees")

def getMiscDir(dataset_dir, family):
  return os.path.join(dataset_dir, "families", family, "misc")

def getRaxmlTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "raxmlGeneTree.newick")

def getRaxmlMultipleTrees(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "raxmlGeneTrees.newick")

def getPhyldogTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "phyldogGeneTree.newick")

def getTreerecsTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "treerecsGeneTree.newick")

def getNotungTree(dataset_dir, family):
  return os.path.join(getTreesDir(dataset_dir, family), "notungGeneTree.newick")

def getALETree(dataset_dir, family, method):
  return os.path.join(getTreesDir(dataset_dir, family), method + "GeneTree.newick")



def getFamilies(dataset_dir):
  return os.listdir(os.path.join(dataset_dir, "families"))

def init_dataset_dir(dataset_dir):
  for family in getFamilies(dataset_dir):
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
  elif (lower_tree == "notung"):
    return os.path.join(gene_trees_dir, "notungGeneTree.newick")
  elif ("GeneRax" in tree):
    return os.path.join(familydir, "results", tree + ".newick")
  elif ("ALE" in tree):
    return os.path.join(gene_trees_dir, tree + "GeneTree.newick")
  elif (lower_tree == "random"):
    return "__random__";
  else:
    return tree

def get_possible_gene_trees():
  return ["raxml", "raxmls", "true", "treerecs", "notung", "phyldog", "random", "ALE-D(T)L", "GeneRax-D(T)L-[Random, Raxml]"]

def get_mapping_file(datadir):
  return os.path.join(datadir, "mapping.link")

def get_alignment_file(datadir):
  return os.path.join(datadir, "alignment.msa")

def get_raxml_model(datadir):
  return os.path.join(datadir, "raxmlBestModel.txt")

def convertToPhyldogSpeciesTree(speciesTree, phyldogSpeciesTree):
  command = "sed s/)[nH][0123456789]*/)/g " + speciesTree #+ " > " + phyldogSpeciesTree
  print(command.split(" "))
  with open(phyldogSpeciesTree, "w") as output:
    subprocess.check_call(command.split(" "), stdout=output)





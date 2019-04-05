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




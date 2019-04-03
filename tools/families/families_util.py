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






import math
import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
sys.path.insert(0, os.path.join("tools", "trees"))
import rescale_bl
import analyze_tree
import experiments as exp
from ete3 import Tree



class SimParameters():
  def set_datadir(self, datadir):
    self.datadir = datadir
    self.simdir = os.path.join(datadir, "lenard")

  def __init__(self):
    self.tag = "highway"
    self.prefix = "lensim"
    self.species_tree = "/data/morelbt/github/lenard/speciesTree.newick"
    self.family_number = 200
    self.loss_rate = 0.1
    self.dup_rate = 0.0
    self.transfer_rate = 0.1
    self.highway_count = 0
    self.highways = [["10", "20"], ["15", "25"]]
    self.highway_rate = 0.2
    self.origination = "root"
    self.sites = 100
    self.max_leaves = 200
    self.blscale = 0.4
    self.set_datadir(fam.get_datadir("lenard_h" + str(self.highway_rate) + "_bl" + str(self.blscale) ))

def run_lenardsim(p):
  command = []
  command.append(exp.lensim_exec)
  command.append("-o")
  command.append(os.path.join(p.simdir, p.prefix))
  command.append("-i")
  command.append(p.origination)
  command.append("-r")
  command.append(str(p.family_number))
  command.append("-d")
  command.append(str(p.dup_rate))
  command.append("-l")
  command.append(str(p.loss_rate))
  command.append("-t")
  command.append(str(p.transfer_rate))
  command.append(p.species_tree)
  if (p.highway_rate != 0.0 and len(p.highways) > 0):
    highwayslist = []
    for highway in p.highways:
      highwayslist.append(highway[0] + "," + highway[1] + "," + str(p.highway_rate)) 
    if (False):
      sources = []
      dests = []
      for highway in p.highways:
        sources.append(highway[0])
        dests.append(highway[1])
      highwaystr = ":".join(sources) + "," + ":".join(dests) + "," + str(p.highway_rate)
    command.append("--highways")
    for h in highwayslist:
      command.append(h)
  print(" ".join(command))
  subprocess.check_call(command)
  
def extract_families(p):
  datadir = p.datadir
  fam.init_top_directories(datadir)
  species_tree = fam.get_species_tree(datadir)
  shutil.copyfile(p.species_tree, species_tree)
  sim_gene_trees = os.path.join(p.simdir, p.prefix + p.origination + ".genetrees")
  index = 0
  for line in open(sim_gene_trees).readlines():
    tree = Tree(line, format = 1)
    leaves = tree.get_leaf_names()
    if (len(leaves) < 4):
      continue
    if (len(leaves) > p.max_leaves):
      continue
    for node in tree.iter_descendants():
      node.dist = node.dist * p.blscale
    family = "family_" + str(index)
    fam.init_family_directories(datadir, family)
    open(fam.get_true_tree(datadir, family), "w").write(tree.write(format = 5))
    mapping = fam.get_mappings(datadir, family)
    with open(mapping, "w") as writer:
      for leaf in leaves:
        gene = leaf
        species = "_".join(leaf.split("_")[:-1])
        writer.write(species + ":" + gene + "\n")
      
    index += 1

def simulate_alignments(datadir, sites):
  seed = 42
  for family in fam.get_families_list(datadir):
    gene_tree = fam.get_true_tree(datadir, family)
    alignment = fam.get_alignment(datadir, family)
    writer = open(alignment, "w")
    command = []
    command.append(exp.seq_gen_exec)
    command.append("-l")
    command.append(str(sites))
    command.append("-m")
    command.append("GTR")
    command.append("-of")
    command.append(gene_tree)
    command.append("-z")
    command.append(str(seed))
    seed += 1
    FNULL = open(os.devnull, 'w')
    with open(alignment, "w") as writer:
      subprocess.check_call(command, stdout=writer, stderr = FNULL)

def generate(p):
  os.mkdir(p.datadir)
  os.mkdir(p.simdir)
  run_lenardsim(p)
  extract_families(p)
  simulate_alignments(p.datadir, p.sites)
  fam.postprocess_datadir(p.datadir)

if (__name__ == "__main__"): 
  p = SimParameters()
  generate(p)


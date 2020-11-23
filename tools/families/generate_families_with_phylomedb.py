import sys
import os
import shutil
import subprocess
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree
from ete3 import Tree
import concurrent.futures
import title_node_names
from itertools import tee
try:
  from itertools import izip
except:
  izip = zip

species_to_remove = set()
species_to_remove.add("HUMAN")
species_to_remove.add("ARATH")

def wget_gunzip(base_url, outdir, filename):
  url = os.path.join(base_url, filename + ".gz")
  outfile = os.path.join(outdir, filename + ".gz")
  wget = []
  wget.append("wget")
  wget.append(url)
  wget.append("--output-document=" + outfile)
  print(" ".join(wget))
  subprocess.check_call(wget)
  subprocess.check_call(["gunzip", outfile])

def get_raw_dir(name):
  return os.path.join(exp.raw_datasets_root, "phylomedb", name)

def dl(index, name):
  rawdir = get_raw_dir(name)
  os.mkdir(rawdir)
  url = "ftp://phylomedb.org/phylomedb/phylomes/phylome_" + str(index)
  wget_gunzip(url, rawdir, "best_trees.txt")
  wget_gunzip(url, rawdir, "phylome_info.txt")

def extract_tree(datadir, family, tree):
  to_keep = []
  unique_names = set()
  for leaf in tree.get_leaves():
    gene = leaf.name
    prefix = gene.split("_")[0]
    species = gene.split("_")[1]
    
    while ((prefix + "_" + species) in unique_names):
      prefix += "o"
    leaf.name = prefix + "_" + species
    unique_names.add(leaf.name)
    if (not species in species_to_remove):
      to_keep.append(leaf.name)
  tree.prune(to_keep, True)
  leaves = tree.get_leaves()
  if (len(leaves) < 4):
    return 
  fam.init_family_directories(datadir, family)
  tree.write(outfile = fam.get_true_tree(datadir, family)) 
  species_to_genes = {}
  for leaf in leaves:
    gene = leaf.name
    species = gene.split("_")[1]
    if (not species in species_to_genes):
      species_to_genes[species] = []
    species_to_genes[species].append(gene)
  mapping_file = fam.get_mappings(datadir, family)
  fam.write_phyldog_mapping(species_to_genes, mapping_file)

def extract(name):
  rawdir = get_raw_dir(name)
  datadir = os.path.join(exp.families_datasets_root, "pdb_" + name)
  fam.init_top_directories(datadir)
  best_trees = os.path.join(rawdir, "best_trees.txt")
  
  for line in open(best_trees).readlines():
    sp = line.split()
    print(sp[0])
    family = sp[0]
    tree = Tree(sp[3], format = 1)
    extract_tree(datadir, family, tree)
  fam.postprocess_datadir(datadir)

def dl_and_extract(index, name):
  dl(index, name)
  extract(name)  



if (__name__ == "__main__"): 
  if (len(sys.argv) < 3): 
    print("Syntax: python " + os.path.basename(__file__) + " db_index outname output")
    exit(1)
  index = sys.argv[1]
  name = sys.argv[2]
  
  dl_and_extract(index, name)


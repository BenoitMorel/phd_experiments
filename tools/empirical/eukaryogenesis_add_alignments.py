import os
import sys
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "trees"))
import experiments as exp
import fam
from read_tree import read_tree
from ete3 import SeqGroup

def subsample_alignment(src, dest, taxa):
  
  msa = SeqGroup(src)
  new_msa = SeqGroup()
  count = 0
  for entry in msa:
    gene = entry[0]
    species = gene.split("_")[0]
    if (species in taxa):
      new_msa.set_seq(entry[0], entry[1])
      count += 1
  assert(count > 0)
  print(count)
  new_msa.write(outfile = dest)
    

def add_alignments(datadir, alignmentdir, species_tree = None):
  leaves = []
  if (species_tree != None):
    leaves = set(read_tree(species_tree).get_leaf_names())
  print("Read an input species tree, will keep " + str(len(leaves)) + " leaves")
  failures = 0
  for ali in os.listdir(alignmentdir):
    print(ali)
    src = os.path.join(alignmentdir, ali)
    family = ali.split(".")[0]
    suffix = "_tokeep"
    path1 = fam.get_family_path(datadir, family + "_PMSF") + suffix
    path2 = fam.get_family_path(datadir, family + "_DEFAULT") + suffix
    if (os.path.isdir(path1)):
      family = family + "_PMSF" + suffix
      assert(not os.path.isdir(path2))
    else:
      family = family  + "_DEFAULT" + suffix
      if (not os.path.isdir(path2)):
        print("Failure " + family)
        failures += 1
        continue

    dest = fam.get_alignment(datadir, family)
    if (len(leaves) > 0):
      subsample_alignment(src, dest, leaves)
    else:
      shutil.copyfile(src, dest)
  print("Number of failures: " + str(failures))

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " datadir alignmentdir species_tree(can be 0).")
    sys.exit(0)
  datadir = sys.argv[1]
  alignmentdir = sys.argv[2]
  species_tree = sys.argv[3]
  if (species_tree == "0"):
    species_tree = None
  add_alignments(datadir, alignmentdir, species_tree)


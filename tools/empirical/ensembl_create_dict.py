import os
import sys
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import fam
import ete3

def create_dict(datadir):
  leaves = ete3.Tree(fam.get_species_tree(datadir)).get_leaf_names()
  with open(fam.get_species_dict(datadir), "w") as writer:
    for leaf in leaves:
      writer.write(leaf)
      writer.write(":")
      writer.write(leaf.replace("UUU", " "))
      writer.write("\n")

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
  datadir = sys.argv[1]
  create_dict(datadir)


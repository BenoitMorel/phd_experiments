import os
import sys
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/families')
import experiments as exp
import fam
from read_tree import read_trees_list



def extract(hogenomfile, datadir):
  families = set(fam.get_families_list(datadir))
  reader = open(hogenomfile)
  family = None
  count_fams = 0
  count_ml = 0
  count_cons = 0
  while (True):
    line = reader.readline()
    if (not line):
      break
    if (line.startswith("AC")):
      family = line.split()[1].replace("\n", "")
      if (not family in families):
        # skip this family
        family = None
      else:
        count_fams += 1
    if (line.startswith("CC   -!- TRIMED TREE Tree in newick format:")):
      reader.readline()
      if (family != None):
        tree = reader.readline()[21:-1]
        open(fam.get_gene_tree_path(datadir, family, "iqtree", "MFP"), "w").write(tree)
        count_ml += 1
    if (line.startswith("CC   -!- TRIMED TREE Consensus tree in newick format:")):
      reader.readline()
      if (family != None):
        tree = reader.readline()[21:-1]
        open(fam.get_gene_tree_path(datadir, family, "iqtreecons", "MFP"), "w").write(tree)
        count_cons += 1
    if (family != None and count_fams % 1000 == 0):
      print("Total families: " + str(len(families)))
      print("Processed families: " + str(count_fams))
      print("ML trees: " + str(count_ml))
      print("Cons trees: " + str(count_cons))


  print("Total families: " + str(len(families)))
  print("Processed families" + str(count_fams))
  print("ML trees: " + str(count_ml))
  print("Cons trees: " + str(count_cons))

if (__name__ == "__main__"):
  if (len(sys.argv) < 3):
    print("Syntax python " + os.path.basename(__file__) + " hogenomfile datadir")
    sys.exit(1)
  hogenomfile = sys.argv[1]
  datadir = sys.argv[2]

  extract(hogenomfile, datadir)

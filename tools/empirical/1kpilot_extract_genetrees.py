"""
  The downloaded 1kpilot trees (for instance from https://datacommons.cyverse.org/browse/iplant/home/shared/onekp_pilot/PNAS_alignments_trees/gene_trees/filtered/filtered_FNA2AA_trees.tar.gz) are organized as a list of directories: "raxmlboot.x.f25" that each contain a gene tree "RAxML_bipartitions.final"
  This scripts output one directory with all gene trees named: "X.newick"
"""

import os
import sys
import shutil

def convert(raw_data_dir, output_tree_dir):
  os.mkdir(output_tree_dir)
  for f in os.listdir(raw_data_dir):
    family = f.split(".")[1]
    src = os.path.join(raw_data_dir, f, "RAxML_bipartitions.final")
    dest = os.path.join(output_tree_dir, family + ".newick")
    try:
      shutil.copyfile(src, dest)
    except:
      print("Can't treat family " + family)


if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " raw_data_dir output_tree_dir")
    sys.exit(1)
  convert(sys.argv[1], sys.argv[2])


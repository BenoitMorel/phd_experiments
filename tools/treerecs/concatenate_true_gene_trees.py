import os
import sys

if (len(sys.argv) != 4):
  print("syntax: python concatenate_true_trees.py trees_dir alignment_file output_file")
  sys.exit(1)

true_dir = sys.argv[1]
alignments_file = sys.argv[2]
output_file = sys.argv[3]

alignments = open(alignments_file).readlines()[1:]

with open(output_file, "w") as f:
  for ali in alignments:
    base = os.path.basename(ali)[0:-16]
    newick = os.path.join(true_dir, base + "tree.nwk")
    print(newick)
    f.write(open(newick).read() + "\n")

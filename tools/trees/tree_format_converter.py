import sys
import os



def newick_to_nexus(input_newick_path, output_nexus_path):
  with open(output_nexus_path, "w") as writer:
    writer.write("#NEXUS\n")
    writer.write("BEGIN TREES;\n")
    writer.write("  TREE tree1 = " + open(input_newick_path).readlines()[0].replace("\n", "") + "\n")
    writer.write("END TREES;\n")



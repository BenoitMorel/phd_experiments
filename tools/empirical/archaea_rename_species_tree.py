import os
import sys
import ete3

def rename(input_tree, species_info, output_tree):
  species_dict = {}
  for line in open(species_info).readlines():
    split = line.split("\t")
    print(split)
    code = split[0]
    long_name = split[1]
    group = ""
    if (len(split) > 2):
      group = " --- " + split[2].split(" ")[0]
    species_dict[code] = long_name + group 
  tree = ete3.Tree(input_tree, format=1)
  for leaf in tree:
    if (leaf.name in species_dict):
      leaf.name = species_dict[leaf.name].replace("\n", "")
  open(output_tree, "w").write(tree.write())
  print(tree.write())

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python archaea_rename_species_tree input_tree species_info output_tree")
    sys.exit(1)
  input_tree = sys.argv[1]
  species_info = sys.argv[2]
  output_tree = sys.argv[3]
  rename(input_tree, species_info, output_tree)


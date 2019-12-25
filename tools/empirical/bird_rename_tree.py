import os
import sys
import ete3
sys.path.insert(0, 'tools/trees')
import analyze_tree

def get_name(split, input_type):
  if (input_type == "short"):
    return split[1]
  elif (input_type == "normal"):
    return split[0]
  elif (input_type == "long"):
    return split[0] + " --- " + split[2] + " --- " + split[3].replace("\n", "")
  elif (input_type == "shortgroup"):
    return split[2]
  elif (input_type == "longgroup"):
    return split[2] + " - " + split[3]
  else:
    print("unknown input type " + input_type)
    sys.exit(1)

def rename(species_info, input_tree, output_tree, input_type, output_type):
  species_dict = {}
  for line in open(species_info).readlines():
    split = " ".join(line.split()).split() # split but also get rid of multiple spaces
    key = get_name(split, input_type)
    value = get_name(split, output_type)
    species_dict[key] = value
  tree = ete3.Tree(input_tree, format=1)
  for leaf in tree:
    if (leaf.name in species_dict):
      leaf.name = species_dict[leaf.name]
  open(output_tree, "w").write(tree.write())
  #print(tree.write())

if (__name__ == "__main__"):
  if (len(sys.argv) != 6):
    print("Syntax python bird_rename_tree.py species_info input_tree output_tree input_type output_type")
    print("type: short, normal, long, shortgroup, longgroup")
    sys.exit(1)
  species_info = sys.argv[1]
  input_tree = sys.argv[2]
  output_tree = sys.argv[3]
  input_type = sys.argv[4]
  output_type = sys.argv[5]
  rename(species_info, input_tree, output_tree, input_type, output_type)
  analyze_tree.count_monophylies(output_tree)


import sys
import os
from read_tree import read_tree


def get_dict(species_tree_path, dictionary_path):
  d = {}
  tree = read_tree(species_tree_path)
  for leaf in tree.get_leaves():
    d[leaf.name] = leaf.name
  if (dictionary_path == None):
    return d
  for line in open(dictionary_path).readlines():
    sp = line.replace("\n", "").split(":")
    d[sp[0]] = sp[1]
  return  d
  
def rename_rec(tree, leaf_dict):
  children = tree.get_children()
  if (len(children) == 0):
    # leaf case
    return leaf_dict[tree.name].split("|")[1:]
  else:
    #internal node
    splits = []
    min_len = 99999
    for child in children:
      sp = rename_rec(child, leaf_dict)
      min_len = min(min_len, len(sp))
      splits.append(sp)
    common_splits = []
    for i in range(0, min_len):
      add = True
      for split in splits[1:]:
        if (split[i] != splits[0][i]):
          add = False
          break
      if (add):
        common_splits.append(splits[0][i])
      else:
        break
  if (len(common_splits) > 0):
    tree.name = common_splits[-1]
  else:
    tree.name = "root"
  return common_splits

def deduplicate(tree, d = {}):
  if (tree.is_leaf()):
    return
  if (not tree.name in d):
    d[tree.name] = 1
  else:
    d[tree.name] += 1
    tree.name = tree.name + "-" + str(d[tree.name])
  for child in tree.get_children():
    deduplicate(child, d)


def rename(input_tree, output_tree, dict_path):
  tree = read_tree(input_tree)
  d = get_dict(input_tree, dict_path)
  rename_rec(tree, d)
  deduplicate(tree)
  tree.write(format = 1, outfile = output_tree)

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__)  + "input_tree output_tree (dict)")
    sys.exit(1)
  input_tree = sys.argv[1]
  output_tree = input_tree
  if (len(sys.argv) > 2):
    output_tree = sys.argv[2]
  dict_path = os.path.join(os.path.dirname(input_tree), os.pardir, "misc", "species_dict.txt")
  if (len(sys.argv) > 3):
    dict_path = sys.argv[3]
  rename(input_tree, output_tree, dict_path)



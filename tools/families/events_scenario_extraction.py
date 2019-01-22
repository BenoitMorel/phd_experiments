import os
import sys
from ete3 import Tree


def read_tree(tree_filename):
  lines = open(tree_filename).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None


def extract_events_from_ale(dataset_dir):

  trees_dir = os.path.join(dataset_dir, "ale", "gene_trees")
  families_dir = os.path.join(dataset_dir, "families")

  events_type = ["S", "D", "T"]
  for tree_filename in os.listdir(trees_dir):
    try:
      tree = read_tree(os.path.join(trees_dir, tree_filename))
    except:
      print("Could not read tree " + os.path.join(trees_dir, tree_filename))
      continue
    family_index = tree_filename.split("_")[1].split(".")[0]
    family_name = family_index + "_pruned"
    events_file = os.path.join(families_dir, family_name, "trueEvents.txt")
    events_count = {}
    for event_type in events_type:
      events_count[event_type] = 0
    for node in tree.traverse("postorder"):
      if (not node.is_leaf()):
        if (not node.name in events_count):
          print("error: unknown event " + node.name)
          exit(1)
        events_count[node.name] += 1
    try:
      writer = open(events_file, "w")
      for event_type in events_type:
        writer.write(event_type + ":" + str(events_count[event_type]) + "\n")
    except:
      print("Ignoring invalid family " + family_name)
      continue

    


if (__name__ == "__main__"):
  dataset_dir = sys.argv[1]
  extract_events_from_ale(dataset_dir)




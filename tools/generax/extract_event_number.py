import sys
import os
import glob
from ete3 import Tree

def get_species_leaves(generaxdir):
  species_tree_file = os.path.join(generaxdir, "species_trees", "starting_species_tree.newick")
  tree = Tree(species_tree_file, format = 1)
  return set(tree.get_leaf_names())

def extract(generaxdir):
  events = {}
  events["S"] = 0
  events["SL"] = 0
  events["D"] = 0
  events["T"] = 0
  events["TL"] = 0
  events["Leaf"] = 0
  f = os.path.join(generaxdir, "per_species_event_counts.txt")
  
  species = get_species_leaves(generaxdir)
  lines = open(f).readlines()
  for line in lines:
    if (line[0] == "#"):
      continue
    sp = line.replace("\n", "").split(" ")
    is_leaf = sp[0] in species
    for event in sp[1:]:
      sp2 = event.split("=")
      key = sp2[0]
      if (key == "S" and is_leaf):
        key = "Leaf"
      events[key] += int(sp2[1])
  return events

def run(generaxdir):
  events = extract(generaxdir)
  for event in events:
    print(event + ":" + str(events[event]))


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " generaxdir")
    sys.exit(1)
  run(sys.argv[1])

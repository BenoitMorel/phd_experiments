import os
import sys
import ete3
import pickle
import subprocess
def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return ete3.Tree(line, format=1)
  return None

def convert_to_notung_tree(input_tree, input_species_tree, mapping_file, notung_tree, notung_species_tree, notung_mapping):
  mapping = {}
  back_mapping = {}
  for line in open(mapping_file).readlines():
    line_split = line.split(":")
    species = line_split[0]
    #mapping[species] = line_split[1][:-1].split(";")
    for gene in line_split[1][:-1].split(";"):
      mapping[gene] = species

  tree = read_tree(input_tree)
  for leaf in tree.get_leaves():
    back_mapping[leaf.name + "_S" + mapping[leaf.name]] = leaf.name
    leaf.name += "_S" + mapping[leaf.name]

  species_tree = read_tree(input_species_tree)
  for leaf in species_tree.get_leaves():
    leaf.name = "S" + leaf.name

  tree.write(format=1, outfile=notung_tree)
  species_tree.write(format=1, outfile=notung_species_tree)
  with open(notung_mapping, 'wb') as f:
    pickle.dump(back_mapping, f, pickle.HIGHEST_PROTOCOL)


def back_convert_notung_tree(notung_tree, notung_mapping, output_tree):
  command = []
  command.append("sed")
  command.append("-e")
  command.append("s/\\[[^][]*\\]//g")
  command.append(notung_tree)
  
  with open(output_tree, "w") as writer:
    subprocess.check_call(command, stdout=writer)
  back_mapping = {}
  with open(notung_mapping, 'rb') as f:
    back_mapping = pickle.load(f)
  tree = read_tree(output_tree)
  good_leaves = []
  for leaf in tree.get_leaves():
    if (not "LOST" in leaf.name):
      good_leaves.append(leaf.name)
  tree.prune(good_leaves)
  for leaf in tree.get_leaves():
    leaf.name = back_mapping[leaf.name]
  tree.write(format=1, outfile=output_tree)
  

if (__name__ == "__main__"):
  if (len(sys.argv) != 7):
    print("Syntax: input_tree input_species_tree  mapping_file notung_tree notung_species_tree notung_mapping")
    exit(1)

  input_tree = sys.argv[1]
  input_species_tree = sys.argv[2]
  mapping_file = sys.argv[3]
  notung_tree = sys.argv[4]
  notung_species_tree = sys.argv[5]
  notung_mapping = sys.argv[6]
  convert_to_notung_tree(input_tree, input_species_tree, mapping_file, notung_tree, notung_species_tree, notung_mapping)
  back_convert_notung_tree("raxmlGeneTree.notung.rearrange.0", notung_mapping, "notungGeneTree.newick")



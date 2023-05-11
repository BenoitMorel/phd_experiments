import os
import sys
import argparse
from ete3 import Tree
import random

palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def parse_arguments(args):
  parser = argparse.ArgumentParser(args)
  parser.add_argument("-i", "--input",
      dest="input",
      help="Input AleRax directory")
  parser.add_argument("-o", "--output",
      dest="output",
      help="Input AleRax directory")
  parser.add_argument("-d", "--label-dictionary",
      dest="species_dict",
      help="Dictionary file to rename the labels")
  op = parser.parse_args()
  return op


def generate_highway(input_path, output_path, d):
  with open(output_path, "w") as writer:
    writer.write("!genetransfer_style\t1\n")
    lines = open(input_path).readlines()
    max_weight = 0.0
    dest_to_lines = {}
# first pass
    for line in lines:
      sp = line.split(",")
      max_weight = max(max_weight, float(sp[0]))
      dest = sp[4].replace("\n", "")
      if (dest in dest_to_lines):
        dest_to_lines[dest].append(line)
      else:
        dest_to_lines[dest] = [line]
       
#second pass
    index = 0
    for highway_dest in dest_to_lines:
      writer.write("# Highways to " + highway_dest + "\n")
      color = palette[index]
      index += 1
      for line in dest_to_lines[highway_dest]:
        sp = line.split(",")
        weight = float(sp[0])
        if (weight < 0.02):
          continue
        support = float(sp[2])
        src = sp[3].replace(" ","")
        dest = sp[4].replace("\n", "")
        src = d[src]
        dest = d[dest]
        writer.write(src + ":" + dest)
        writer.write("\tdir=0:1,color=" + color + ",linewidth=" + str(5.0 * weight / max_weight) + "\n")


def get_label_dict(tree, dictionary_path):
  d = {}
  if (dictionary_path == ""):
    dictionary_path = None
  
  for leaf in tree.get_leaves():
    d[leaf.name] = leaf.name
  if (None == dictionary_path):
    return d
  if (not os.path.isfile(dictionary_path)):
    print("Label dictionary file " + dictionary_path + " not found")
    sys.exit(1)
  for line in open(dictionary_path).readlines():
    sp = line.replace("\n", "").split(":")
    d[sp[0]] = sp[1]
  return  d

def get_species_dict(species_tree, species_dict = None):
  tree = Tree(species_tree, format = 1)
  label_dict = get_label_dict(tree, species_dict)
  d = {}
  d_left = {}
  d_right = {}
  for node in tree.traverse("postorder"):
    if (node.is_leaf()):
      d[node.name] = label_dict[node.name]
      node.name = d[node.name]
      d_left[node.name] = node.name
      d_right[node.name] = node.name
    else:
      children = node.children
      if (len(children) != 2):
        print("Error, the species tree is not binary")
        sys.exit(1)
      left = children[0]
      right = children[1]
      # evolview internal node naming
      d[node.name] = d_left[left.name] + "," + d_right[right.name]
      d_left[node.name] = d_left[left.name]
      d_right[node.name] = d_right[right.name]
  return d

def generate(op):
  try:
    os.mkdir(op.output)
  except:
    print("Can't create directory " + op.input + ". Please check that the directory does not exist yet")
    sys.exit(1)
  
  print("Initialization...")
  species_tree = os.path.join(op.input, "species_trees", "inferred_species_tree.newick")
  d = get_species_dict(species_tree, op.species_dict)
  highways_path = os.path.join(op.input, "highway_accepted_highways.txt")
  highways_evolview = os.path.join(op.output, "highways.txt")
  if (os.path.isfile(highways_path)):
    print("Highways file detected, generating highways evolview file " + highways_evolview)
    generate_highway(highways_path, highways_evolview, d)


if (__name__== "__main__"):
  op = None
  try:
    op = parse_arguments(sys.argv)
  except Exception as inst:
    print("[Error] " + str(type(inst)) + " " + str(inst)) 
    sys.exit(1)
  generate(op)

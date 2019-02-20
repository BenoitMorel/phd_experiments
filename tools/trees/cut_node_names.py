import re
import sys

def shorten_obj_keep_first(o):
  return o.group(0).split("_")[0]

def shorten_obj_remove_last(o):
  return "_".join(o.group(0).split("_")[:-1])

def cut_node_names(input_tree, output_tree, keep_first):
  s = open(input_tree).read()
  if (keep_first):
    res = re.sub("[A-Za-z0-9_]*", shorten_obj_keep_first, s)
  else:
    res = re.sub("[A-Za-z0-9_]*", shorten_obj_remove_last, s)
  with open(output_tree, "w") as writer:
    writer.write(res)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: input_tree output_tree")
    exit(1)
  cut_node_names(sys.argv[1], sys.argv[2])





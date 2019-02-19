import re
import sys

def shorten_obj(o):
  return o.group(0).split("_")[0]

def cut_node_names(input_tree, output_tree):
  s = open(input_tree).read()
  res = re.sub("[A-Za-z0-9_]*", shorten_obj, s)
  print(s)
  print(res)
  with open(output_tree, "w") as writer:
    writer.write(res)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: input_tree output_tree")
    exit(1)
  cut_node_names(sys.argv[1], sys.argv[2])





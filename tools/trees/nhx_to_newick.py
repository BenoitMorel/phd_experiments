import re
import sys

def nhx_to_newick(input_tree, output_tree):
  s = open(input_tree).read()
  res = re.sub("\[.*?\]", "", s) + ";"
  with open(output_tree, "w") as writer:
    writer.write(res)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: input_tree output_tree")
    exit(1)
  nhx_to_newick(sys.argv[1], sys.argv[2])




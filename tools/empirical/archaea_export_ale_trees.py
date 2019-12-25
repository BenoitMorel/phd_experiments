import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import experiments as exp

def get_tree(ale_file):
  for line in open(ale_file):
    if (";" in line):
      return line
  return None

def export_ale_trees(input_trees_dir, datadir):
  for f in os.listdir(input_trees_dir):
    family = f.replace("fix_", "").replace("_c60_lg_combined.ale", "")
    family = family + "_bmge30_renamed"
    output_tree = fam.build_gene_tree_path(datadir, "c60", family, "ale")
    ale_file = os.path.join(input_trees_dir, f)
    tree_string = get_tree(ale_file)
    ok = True
    if (tree_string):
      print("Not tree for " + family)
    try:
      with open(output_tree, "w") as writer:
        writer.write(tree_string)
    except:
      ok = False
    if (not ok):
      print("Error with family " + family)

      

if (__name__ == "__main__"): 
  if (len(sys.argv) != 3):
    print("syntax: input_trees_dir datadir")
    exit(1)

  input_trees_dir = sys.argv[1]
  datadir = sys.argv[2]
  export_ale_trees(input_trees_dir, datadir)



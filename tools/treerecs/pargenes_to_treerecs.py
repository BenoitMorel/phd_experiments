import os
import sys
from subprocess import call

def get_alignment_dir(pargenes_dir):
  pargenes_logs_file = os.path.join(pargenes_dir, "pargenes_logs.txt")
  logs_lines = open(pargenes_logs_file).readlines()
  pargenes_arguments = []
  alignment_dir = ""
  for i in range(0, len(logs_lines)):
    if "ParGenes was called as follow" in logs_lines[i]:
      pargenes_arguments = logs_lines[i+1].split(" ")
      break
  for i in range(0, len(pargenes_arguments)):
    if pargenes_arguments[i] in ["-a", "--alignments-dir"]:
      alignment_dir = pargenes_arguments[i + 1]
      break
  return alignment_dir



if (len(sys.argv) != 3):
  print("Syntax: python pargenes_to_treerecs.py pargenes_dir treerecs_dir")
  sys.exit(1)
pargenes_dir = sys.argv[1]
treerecs_dir = sys.argv[2]

if (os.path.isdir(treerecs_dir)):
  print("Error: directory " + treerecs_dir + " already exists")
  sys.exit(1)

os.makedirs(treerecs_dir)
alignment_dir = get_alignment_dir(pargenes_dir)
print("alignment dir " + alignment_dir)
msa_list = os.listdir(alignment_dir)

alignments_path = os.path.join(treerecs_dir, "alignment.txt")
gene_trees_path = os.path.join(treerecs_dir, "geneTrees.newick")

print(alignments_path)

alignments = open(alignments_path, "w")
gene_trees = open(gene_trees_path, "w")

alignments.write("GTR\n")

support_dir = os.path.join(pargenes_dir, "supports_run", "results")



for msa in msa_list:
  msa_name = msa.replace(".", "_")
  support_path = os.path.join(support_dir, msa_name + ".support.raxml.support")
  if not os.path.isfile(support_path):
    print("no support file for " + msa)
    continue
  alignments.write(os.path.join(alignment_dir, msa) + "\n")
  gene_trees.write(open(os.path.join(support_path)).read())

alignments.close()
gene_trees.close()

script_dir = os.path.dirname(os.path.realpath(__file__))
call(["cp", gene_trees_path, gene_trees_path + "_old"])
command = [os.path.join(script_dir, "raxml_to_treerecs_supporttree.sh"), gene_trees_path]
print("command " + str(command))
call(command)



import os
import sys




if (len(sys.argv) != 5):
  print("Syntax: python alignment_splitter.py input_alignment input_gene_trees output_dir prefix_path")
  sys.exit(1)

input_alignment = sys.argv[1]
input_gene_trees = sys.argv[2]
output_dir = sys.argv[3]
prefix_path = sys.argv[4]

if (os.path.exists(output_dir)):
  print("output dir already exists")
  sys.exit(1)

os.makedirs(output_dir)
ali_dir = os.path.join(output_dir, "split_alignments")
trees_dir = os.path.join(output_dir, "split_gene_trees")
os.makedirs(ali_dir)
os.makedirs(trees_dir)


trees = open(input_gene_trees).readlines()

with open(input_alignment) as f:
  model = f.readline()
  i = 0
  for line in f:

    line = line[:-1]
    base = os.path.basename(line).split(".")[0]
    with open(os.path.join(ali_dir, base + ".ali"), "w") as writer:
        writer.write(model)
        writer.write(os.path.join(prefix_path, line))
    with open(os.path.join(trees_dir, base + ".newick"), "w") as writer:
      writer.write(trees[i])
    i += 1


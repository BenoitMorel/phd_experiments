import os
import sys




if (len(sys.argv) != 4):
  print("Syntax: python alignment_splitter.py input_alignment output_dir prefix_path")
  sys.exit(1)

input_alignment = sys.argv[1]
output_dir = sys.argv[2]
prefix_path = sys.argv[3]

if (os.path.exists(output_dir)):
  print("output dir already exists")
  sys.exit(1)

os.makedirs(output_dir)

with open(input_alignment) as f:
  model = f.readline()
  for line in f:
    line = line[:-1]
    base = os.path.basename(line).split(".")[0]
    with open(os.path.join(output_dir, base + ".ali"), "w") as writer:
        writer.write(model)
        writer.write(os.path.join(prefix_path, line))



import os
import sys

def generate(tsv_file, output):
  lines = open(tsv_file).readlines()[1:]
  with open(output, "w") as writer:
    for line in lines:
      sp = line.split()
      print(sp)
      taxa = sp[0 ]
      sp2 = sp[2].split(";")
      long_name = taxa + "_" + sp2[1] + "_" + sp2[2]
      writer.write(taxa + ":")
      writer.write(long_name)
      writer.write("\n")

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python " + os.path.basename(__file__) + " tsv_file output")
    sys.exit(1)
  tsv_file = sys.argv[1]
  output = sys.argv[2]
  generate(tsv_file, output)


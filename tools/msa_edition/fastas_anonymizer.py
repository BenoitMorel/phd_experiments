import sys
import os

if (len(sys.argv) != 2):
  print("Syntax: python fastas_anonymizer.py inputdir outputdir")

inputdir = sys.argv[1]
outputdir = sys.argv[2]
if (os.path.isdir(outputdir)):
  print("error: output directory already exsits")
  sys.exit(1)

os.makedirs(outputdir)

input_files = os.listdir(inputdir)

for fasta in input_files:
  print("treating " + fasta)
  with open(os.path.join(outputdir, fasta), "w") as w:
    lines = open(os.path.join(inputdir, fasta)).readlines() 
    index = 0
    for line in lines:
      if (line.startswith(">")):
        line = ">seq" + str(index) + "\n"
        index += 1
      w.write(line)


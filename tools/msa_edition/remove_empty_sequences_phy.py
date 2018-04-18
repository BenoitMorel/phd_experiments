import os
import sys


def remove_empty_sequences(inputFile, outputFile):
  lines = open(inputFile, "r").readlines()
  sites = int(lines[0].split(" ")[1][:-1])
  empty = "-" * sites
  filtered_lines = []
  for line in lines:
    seq = line.split(" ")[1][:-1]
    if (seq != empty):
      filtered_lines.append(line)
  filtered_lines[0] = str(len(filtered_lines) - 1) + " " + str(sites) + "\n"
  with open(outputFile, "w") as writer:
    for line in filtered_lines:
      writer.write(line)

inputDir = sys.argv[1]
outputDir = sys.argv[2]

if (inputDir == outputDir):
  print("Error: input and output directories are the same")
  sys.exit(1)

try:
  os.makedirs(outputDir)
except:
  pass

for base in os.listdir(inputDir):
  inputFile = os.path.join(inputDir, base)
  outputFile = os.path.join(outputDir, base)
  remove_empty_sequences(inputFile, outputFile)

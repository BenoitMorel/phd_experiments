import os
import sys
import re
import ete3

def fullmatch(regex, string, flags=0):
  return None != re.match("(?:" + regex + r")\Z", string, flags=flags)

def is_empty_dna(sequence):
  pattern = "[^AGTCagct]*"
  return fullmatch(pattern, sequence)

"""
  return a new MSA with all empty sequences removed
  An empty sequence is a sequence that does not countain
  any ACGTacgt
"""
def get_cleaned_msa_dna(msa):
  new_msa = ete3.SeqGroup()
  for entry in msa.get_entries():
    if (not is_empty_dna(entry[1])):
      new_msa.set_seq(entry[0], entry[1])
  return new_msa

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
if (__name__ == "__main__"): 

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

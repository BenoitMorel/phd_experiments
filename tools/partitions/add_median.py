import sys
import os

"""
Open a partition file, and add an A after each G[...] 
to select median gamma rates instead of mean.
Then output the new partitions into another file
"""
if (len(sys.argv) != 3):
  print("usage: python add_median.py input_part output_part")
  sys.exit(1)
input_part = sys.argv[1]
output_part = sys.argv[2]

if (input_part == output_part):
  print("Error: both files are the same. Exiting")
  sys.exit(2)

lines = open(input_part).readlines()
with open(output_part, "w") as writer:
  for line in lines:
    split_line = line.split(",")
    
    split = split_line[0].split("+")
    for idx, elem in enumerate(split):
      if (elem.startswith("G")):
        split[idx] = elem + "A"
    split_line[0] = "+".join(split)
    writer.write(",".join(split_line))


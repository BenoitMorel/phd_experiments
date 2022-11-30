import sys

def convert(treerecs_filename, phyldog_filename):
  ilines = open(treerecs_filename).readlines()
  d = {}
  for line in ilines:
    split = line[:-1].split()
    gene = split[0]
    species = split[1]
    if (not species in d):
      d[species] = []
    d[species].append(gene)
  with open(phyldog_filename, "w") as writer:
    for species in d:
      writer.write(species + ":" + ";".join(d[species]) + "\n")


if (__name__ == "__main__"):

  if (3 != len(sys.argv)):
    print("Syntax: python " + os.path.basename(__file__) + "treerecs_filename phyldog_filename")
    sys.exit()
  convert(sys.argv[1], sys.argv[2])



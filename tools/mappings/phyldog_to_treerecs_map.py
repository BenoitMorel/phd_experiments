import sys

def convert(phyldog_filename, treerecs_filename):
  ilines = open(phyldog_filename).readlines()
  with open(treerecs_filename, "w") as f:
    for line in ilines:
      split = line[:-1].split(":")
      species = split[0]
      genes = split[1].split(";")
      for gene in genes:
        f.write(gene + " " + species + "\n")

if (__name__ == "__main__"):

  if (3 != len(sys.argv)):
    print("Syntax: python " + os.path.basename(__file__) + " phyldog_filename treerecs_filename")
    sys.exit()
  convert(sys.argv[1], sys.argv[2])


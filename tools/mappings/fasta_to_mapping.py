import sys

def fasta_to_mapping(input_fasta, output_mapping):
  lines = open(input_fasta).readlines()
  dico = {}
  for line in lines:
    if (not line.startswith(">")):
      continue
    split = line[1:].split("_")
    gene = split[0]
    species = split[1][:-1]
    if (not species in dico):
      dico[species] = []
    dico[species].append(gene + "_" + species)

  with open(output_mapping, "w") as writer:
    for species in dico:
      writer.write(species + ":" + ";".join(dico[species]) + "\n")



if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: python phy_to_fasta.py input_fasta output_mapping")
    exit(1)
  fasta_to_mapping(sys.argv[1], sys.argv[2])


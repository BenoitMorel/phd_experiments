import os
import sys

sys.path.insert(0, 'tools/msa_edition')
import remove_empty_sequences

is_dna = True

def extract_msa(s, curs, max_curs, writer):
  while (curs < max_curs and curs != 0):
    end_species = s.find("\t", curs)
    species = s[curs: end_species]
    curs_seq = end_species + 1
    end_seq = s.find("\n", curs_seq)
    seq = s[curs_seq:end_seq]
    curs = s.find("\t", end_seq)
    curs = curs + 1
    if (remove_empty_sequences.is_empty(seq, is_dna)):
      continue
    writer.write(">")
    writer.write(species)
    writer.write("\n")
    writer.write(seq)
    writer.write("\n")

def extract_msas(nexus, outputdir):
  os.mkdir(outputdir)
  s = open(nexus).read()
  curs = 0
  count = 0
  while (True):
    # search for family name
    curs = s.find(":", curs)
    if (curs == -1):
      break
    end_family = s.find("]", curs)
    family = s[curs + 2 : end_family]
    fasta = os.path.join(outputdir, family + ".fasta")
    print(fasta)
    curs = s.find("[", end_family)
    if (curs == -1):
      curs = len(s)
    with open(fasta, "w") as writer:
      extract_msa(s, end_family + 4, curs, writer) 
    count = count + 1

  print(count)


if (__name__ == "__main__"):
  if (len(sys.argv) < 3):
    print("Syntax python " + os.path.basename(__file__) + " nexus outputdir")
    sys.exit(1)
  
  print("WARNING: the script assumes DNA data")
  nexus = sys.argv[1]
  outputdir = sys.argv[2]
  extract_msas(nexus, outputdir) 
  print("WARNING: the script assumes DNA data")




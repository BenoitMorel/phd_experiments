import os
import sys

class Gene:
  def __init__(self):
    self.species = ""
    self.gene = ""
    self.chrom = ""
    self.begin = 0
    self.end = 0
    self.family = ""
    self.before = None
    self.after = None

class EmfInfo:
  def __init__(self):
    chromosomes = {}
    allgenes = {}
    families = {}


def add_to_chromosomes(gene, chromosomes):
  if (not gene.species in chromosomes):
    chromosomes[gene.species] = {}
  if (not gene.chrom in chromosomes[gene.species]):
    chromosomes[gene.species][gene.chrom] = []
  chromosomes[gene.species][gene.chrom].append(gene)




def read(emf):
  families = {}
  chromosomes = {}
  allgenes = {}

  index = 1
  family = "fam" + str(index)
  families[family] = []

  for line in open(emf).readlines():
    if (line.startswith("DATA")):
      index += 1
      family = "fam" + str(index)
      families[family] = []
    if (not line.startswith("SEQ")):
      continue
    sp = line.split()
    g = Gene()
    g.species = sp[1]
    g.gene = sp[2]
    g.chrom = sp[3]  
    g.begin = int(sp[4])
    g.end = int(sp[5])
    g.family = family
    families[family].append(g)
    add_to_chromosomes(g, chromosomes)
    allgenes[g.gene] = g
   
  for species in chromosomes:
    for chrom_id in chromosomes[species]:
      chromosome = chromosomes[species][chrom_id]
      chromosome.sort(key = lambda x: x.begin)
      for i in range(0, len(chromosome) - 1):
        gene1 = chromosome[i]
        gene2 = chromosome[i+1]
        gene1.after = gene2.gene
        gene2.before = gene1.gene


  info = EmfInfo()
  info.families = families
  info.chromosomes = chromosomes
  info.allgenes = allgenes
  return info


if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + " emf output_file")
  emf = sys.argv[1]
  read(emf)


import sys
import os
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree
from ete3 import Tree
import concurrent.futures

class SeqEntry():
  def __init__(self, line, family):
    split = line.split(" ")
    temp = split[1].split("_")
    self.species = temp[0].title() + "_" + temp[1].title()
    self.gene = split[2]
    self.chromozome = split[3]
    self.begin = int(split[4])
    self.end = int(split[5])
    self.wtf1 = split[6]
    self.wtf2 = split[7]
    self.wtf_function = split[8][:-1]
    self.family = family

  def __str__(self):
    res = "("
    res += "species:" + self.species + ", "
    res += "gene:" + self.gene + ", "
    res += "chrom:" + self.chromozome + ", "
    res += "begin:" + str(self.begin) + ", "
    res += "end:" + str(self.end) + ", "
    res += "family:" + self.family
    res += ")"
    return res

def parse_nhx_emf(emf_file, species_dict):
  seq_entries_dict = {}  
  family_index = 0
  family = "family_" + str(family_index)
  for line in open(emf_file).readlines():
    if (line[:3] == "SEQ"):
      entry = SeqEntry(line, family)
      if (entry.species in species_dict):
        seq_entries_dict[entry.gene] = entry
    elif (line[:4] == "DATA"):
      family_index += 1
      family = "family_" + str(family_index)
  print("number of seq entries: " + str(len(seq_entries_dict)))
  return seq_entries_dict


def parse_fasta(fasta_file, genes_dict):
  alignments_dico = {}
  cur_name = ""
  cur_seq = ""
  for line in open(fasta_file).readlines():
    if (line.strip() == "" or line[0] == "/"):
      continue
    if (line[0] == ">"):
      if (len(cur_name) and cur_name in genes_dict):
        alignments_dico[cur_name] = cur_seq
      cur_name = line[1:-1]
      cur_seq = line
    else:
      cur_seq += line
  alignments_dico[cur_name] = cur_seq
  print("number of DNA seq: " + str(len(alignments_dico)))
  return alignments_dico

def parallelized_function(params):
  datadir, family = params
  create_random_tree.create_random_tree(fam.getAlignment(datadir, family), fam.getTrueTree(datadir, family))

# this function should be common to all generation scripts
def prepare_datadir(datadir):
  # phyldog species trees
  fam.convertToPhyldogSpeciesTree(fam.get_species_tree(datadir), fam.get_phyldog_species_tree(datadir)) 
  exp.mkdir(os.path.join(datadir, "alignments"))
  # alignments
  for family in fam.getFamiliesList(datadir):
    family_dir = fam.getFamily(datadir, family)
    exp.relative_symlink(fam.getAlignment(datadir, family), os.path.join(datadir, "alignments", family + ".fasta"))
    fam.convert_phyldog_to_treerecs_mapping(fam.get_mappings(datadir, family), fam.get_treerecs_mappings(datadir, family)) 
    exp.relative_symlink(fam.getSpeciesTree(datadir), os.path.join(fam.getFamiliesDir(datadir), family, "speciesTree.newick"))
  print("create random trees")
  params = []
  for family in fam.getFamiliesList(datadir):
    params.append((datadir, family))
  with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(parallelized_function, params)

def export_msa(seq_entries, alignments_dico, output_file):
  with open(output_file, "w") as writer:
    for seq_entry in seq_entries:
      writer.write(alignments_dico[seq_entry.gene])

def export_mappings(families_seq_entries, output_file):
  species_to_genes = {}
  for entry in families_seq_entries:
    if (not entry.species in species_to_genes):
      species_to_genes[entry.species] = []
    species_to_genes[entry.species].append(entry.gene)
  fam.write_phyldog_mapping(species_to_genes, output_file)

def export(species_tree, seq_entries_dict, alignments_dico, datadir):
  per_family_seq_entries = {}
  for gene in seq_entries_dict:
    seq_entry = seq_entries_dict[gene]
    if (not seq_entry.family in per_family_seq_entries):
      per_family_seq_entries[seq_entry.family] = []
    per_family_seq_entries[seq_entry.family].append(seq_entry)
  
  os.makedirs(datadir)
  shutil.copy(species_tree, fam.getSpeciesTree(datadir))
  
  families_dir = fam.getFamiliesDir(datadir)
  os.makedirs(families_dir)
  print("Number of families: " + str(len(per_family_seq_entries)))
  for family in per_family_seq_entries:
    family_dir = fam.getFamily(datadir, family)
    os.makedirs(family_dir)
    seq_entries = per_family_seq_entries[family]
    export_msa(seq_entries, alignments_dico, fam.getAlignment(datadir, family)) 
    export_mappings(seq_entries, fam.get_mappings(datadir, family))
  prepare_datadir(datadir)

def get_species_dict(species_tree_file):
  res = {}
  tree = Tree(species_tree_file, 1)
  for leaf in tree.get_leaves():
    split = leaf.name.split("_")
    res[split[0].title() + "_" + split[1].title()] = True
  return res

def extract_from_ensembl(nhx_emf_file, fasta_file, species_tree, datadir):
  species_dict = get_species_dict(species_tree)
  seq_entries_dict = parse_nhx_emf(nhx_emf_file, species_dict)
  alignments_dico = parse_fasta(fasta_file, seq_entries_dict)
  export(species_tree, seq_entries_dict, alignments_dico, datadir)

if (__name__ == "__main__"): 
  if (len(sys.argv) != 5): 
    print("Syntax: python " + os.path.basename(__file__) + " nhx_emf fasta output")
    exit(1)
  nhx_emf_file = sys.argv[1]
  fasta_file = sys.argv[2]
  species_tree = sys.argv[3]
  datadir = sys.argv[4]
  extract_from_ensembl(nhx_emf_file, fasta_file, species_tree, datadir)
  

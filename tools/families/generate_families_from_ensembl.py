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
import title_node_names
from itertools import tee, izip

class SeqEntry():
  def __init__(self, line, family):
    split = line.split(" ")
    self.species = title_node_names.get_title(split[1])
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

def compare(entry1, entry2):
  if (entry1.species != entry2.species):
    return cmp(entry1.species, entry2.species)
  if (entry1.chromozome != entry2.chromozome):
    return cmp(entry1.chromozome, entry2.chromozome)
  return cmp(entry1.begin, entry2.begin)

# read an ensembl tree, prune the genes that we do not consider,
# and return true if the resulting tree has more than 3 taxa
def read_ensembl_tree(line, family, trees_dict, seq_entries_dict):
  tree = Tree(line.replace("\n", ""), format = 1)
  leaves = tree.get_leaves()
  filtered_leaves = []
  for leaf in leaves:
    if leaf.name in seq_entries_dict:
      filtered_leaves.append(leaf.name)
  try:
    tree.prune(filtered_leaves)
  except:
    return False
  print(len(tree.get_leaves()))
  if (len(tree.get_leaves()) > 3):
    trees_dict[family] = line.replace("\n", "")
    return True
  return False


def parse_nhx_emf(emf_file, species_dict, max_families):
  seq_entries_dict = {}  
  current_entries_dict = {}
  trees_dict = {}
  family_index = 0
  family = "family_" + str(family_index)
  last_was_data = False
  for line in open(emf_file).readlines():
    if (last_was_data):
      last_was_data = False
      if (read_ensembl_tree(line, family, trees_dict, current_entries_dict)):
        for gene in current_entries_dict:
          seq_entries_dict[gene] = current_entries_dict[gene]
        family_index += 1
        family = "family_" + str(family_index)
        print(family)
      current_entries_dict = {}
      if (family_index >= max_families):
        break
    if (line[:3] == "SEQ"):
      entry = SeqEntry(line, family)
      if (entry.species in species_dict):
        current_entries_dict[entry.gene] = entry
    elif (line[:4] == "DATA"):
      last_was_data = True

  print("number of seq entries: " + str(len(seq_entries_dict)))
  return seq_entries_dict, trees_dict


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
  #print("create random trees")
  #params = []
  #for family in fam.getFamiliesList(datadir):
  #  params.append((datadir, family))
  #with concurrent.futures.ProcessPoolExecutor() as executor:
  #  executor.map(parallelized_function, params)

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

def pairwise(iterable):
  "s -> (s0,s1), (s1,s2), (s2, s3), ..."
  a, b = tee(iterable)
  next(b, None)
  return izip(a, b)


def extract_adjacencies(seq_entries_dict, datadir)):
  seq_entries_list = []
  for gene in seq_entries_dict:
    seq_entries_list.append(seq_entries_dict[gene])
  seq_entries_list.sort(lambda x,y: compare(x, y))
  with open(fam.get_adjacencies(datadir), "w") as writer:
    for gene1, gene2 in pairwise(seq_entries_list):
      if (gene1.species == gene2.species and gene1.chromozome == gene2.chromozome):
        writer.write(gene1.gene + "|" + gene2.gene + "\n")
  with open(fam.get_prefixed_adjacencies(datadir), "w") as writer:
    for gene1, gene2 in pairwise(seq_entries_list):
      if (gene1.species == gene2.species and gene1.chromozome == gene2.chromozome):
        writer.write(gene1.species + "_" + gene1.gene + "|" + gene2.species + "_" + gene2.gene + "\n")

def extract_deco_mappings(seq_entries_dict, output_file):
  with open(output_file, "w") as writer:
    for gene in seq_entries_dict:
      writer.write(seq_entries_dict[gene].species + " " + gene + "\n")

def export(species_tree, seq_entries_dict, trees_dict, alignments_dico, datadir):
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
  processed_families = 0
  for family in per_family_seq_entries:
    family_dir = fam.getFamily(datadir, family)
    os.makedirs(family_dir)
    seq_entries = per_family_seq_entries[family]
    export_msa(seq_entries, alignments_dico, fam.getAlignment(datadir, family)) 
    export_mappings(seq_entries, fam.get_mappings(datadir, family))
    processed_families += 1
    with open(fam.getTrueTree(datadir, family), "w") as writer:
      writer.write(trees_dict[family])
  extract_adjacencies(seq_entries_dict, datadir)
  extract_deco_mappings(seq_entries_dict, fam.get_deco_mappings(datadir))
  prepare_datadir(datadir)

def get_species_dict(species_tree_file):
  res = {}
  tree = Tree(species_tree_file, 1)
  for leaf in tree.get_leaves():
    split = leaf.name.split("_")
    res[leaf.name] = True
  return res


def extract_from_ensembl(nhx_emf_file, fasta_file, species_tree, datadir, max_families):
  new_species_tree = species_tree + ".titled"
  title_node_names.title_node_names(species_tree, new_species_tree)
  species_tree = new_species_tree
  species_dict = get_species_dict(species_tree)
  print(species_dict)
  seq_entries_dict, trees_dict = parse_nhx_emf(nhx_emf_file, species_dict, max_families)
  
  alignments_dico = parse_fasta(fasta_file, seq_entries_dict)
  export(species_tree, seq_entries_dict, trees_dict, alignments_dico, datadir)


def get_max_families(argv):
  max_families = 10000000000
  for i in range(0, len(argv)):
    if (argv[i] == "--max-families"):
      max_families = int(argv[i+1])
      print("reading max families " + str(max_families))
  return max_families

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5): 
    print("Syntax: python " + os.path.basename(__file__) + " nhx_emf fasta output")
    exit(1)
  nhx_emf_file = sys.argv[1]
  fasta_file = sys.argv[2]
  species_tree = sys.argv[3]
  datadir = sys.argv[4]
  max_families = get_max_families(sys.argv)
  extract_from_ensembl(nhx_emf_file, fasta_file, species_tree, datadir, max_families)
  

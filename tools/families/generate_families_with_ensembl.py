import sys
import os
import shutil
import functools
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import experiments as exp
import fam
import create_random_tree
from ete3 import Tree
import concurrent.futures
import title_node_names
from itertools import tee
try:
  from itertools import izip
except:
  izip = zip


species_to_remove = set()
species_to_remove.add("Oryzias_latipes_ASM223467v1")
species_to_remove.add("Oryzias_latipes_ASM223471v1")
species_to_remove.add("Oryzias_latipes_ASM223469v1")
species_to_remove.add("Bos_indicus_x_Bos_taurus_UOA_Brahman_1")
species_to_remove.add("Bos_indicus_x_Bos_taurus_UOA_Angus_1")
species_to_remove.add("Mus_musculus_strain_reference_CL57BL6")
species_to_remove.add("Mus_musculus_strain_129S1/SvImJ")
species_to_remove.add("Mus_musculus_strain_NOD/ShiLtJ")
species_to_remove.add("Mus_musculus_strain_A/J")
species_to_remove.add("Mus_musculus_strain_C3H/HeJ")
species_to_remove.add("Mus_musculus_domesticus")
species_to_remove.add("Mus_musculus_strain_DBA/2J")
species_to_remove.add("Mus_musculus_musculus")
species_to_remove.add("Mus_musculus_strain_NZO/HlLtJ")
species_to_remove.add("Mus_musculus_strain_C57BL/6NJ")
species_to_remove.add("Mus_musculus_strain_LP/J")
species_to_remove.add("Mus_musculus_strain_CBA/J")
species_to_remove.add("Mus_musculus_strain_FVB/NJ")
species_to_remove.add("Mus_musculus_strain_BALB/cJ")
species_to_remove.add("Mus_musculus_strain_AKR/J")
species_to_remove.add("Mus_musculus_castaneus")
#species_to_remove.add("Sus_scrofa_strain_reference")
species_to_remove.add("Sus_scrofa_strain_rongchang")
species_to_remove.add("Sus_scrofa_strain_pietrain")
species_to_remove.add("Sus_scrofa_strain_meishan")
species_to_remove.add("Sus_scrofa_strain_bamei")
species_to_remove.add("Sus_scrofa_strain_largewhite")
species_to_remove.add("Sus_scrofa_strain_berkshire")
species_to_remove.add("Sus_scrofa_strain_wuzhishan")
species_to_remove.add("Sus_scrofa_strain_usmarc")
species_to_remove.add("Sus_scrofa_strain_hampshire")
species_to_remove.add("Sus_scrofa_strain_tibetan")
species_to_remove.add("Sus_scrofa_strain_landrace")
species_to_remove.add("Sus_scrofa_strain_jinhua")



species_to_replace = {}
# VERTEBRATES
species_to_replace["Gorilla_gorilla_gorilla"] = "Gorilla_gorilla"
species_to_replace["Mus_caroli_strain_CAROLI_EIJ"] = "Mus_Caroli"
species_to_replace["Sus_scrofa_strain_reference"] = "Sus_Scrofa"
species_to_replace["Mus_spretus_strain_SPRET/EiJ"] = "Mus_Spretus"
species_to_replace["Mus_pahari_strain_PAHARI_EIJ"] = "Mus_Pahari"
species_to_replace["Cricetulus_griseus"] = "Cricetulus_Griseus_Chok1Gshd"
species_to_replace["Cricetulus_griseus_CHOK1GS_HDv1"] = "Cricetulus_Griseus_Chok1Gshd"
species_to_replace["Cricetulus_griseus_CriGri_1.0"] = "Cricetulus_Griseus_CriGri"
species_to_replace["Cricetulus_griseus_CriGri-PICR"] = "Cricetulus_Griseus_PICR"
species_to_replace["Astyanax_mexicanus_Astyanax_mexicanus-2.0"] = "astyanax_mexicanus"
species_to_replace["Astyanax_mexicanus_Astyanax_mexicanus-1.0.2"] = "astyanax_mexicanus_pachon"
species_to_replace["Canis_lupus_familiaris"] = "Canis_familiaris"
species_to_replace["Heterocephalus_glaber_HetGla_1.0"] = "heterocephalus_glaber_male"
species_to_replace["Heterocephalus_glaber_HetGla_female_1.0"] = "heterocephalus_glaber_female"
species_to_replace["Cebus_capucinus_imitator"] = "Cebus_capucinus"
species_to_replace["Saccharomyces_cerevisiae_strain_S288C"] = "Saccharomyces_cerevisiae"
species_to_replace["Caenorhabditis_elegans_strain_N2"] = "Caenorhabditis_elegans"
#species_to_replace["Mus_musculus_strain_reference_CL57BL6"] = "Mus_Musculus"







def load_species_dict(filepath):
  d = {}
  for line in open(filepath).readlines():
    sp = line.replace("\n", "").split("=")
    if (len(sp) == 2):
      d[sp[0]] = sp[1]
  return d

def load_species_blacklist(filepath):
  b = set()
  for line in open(filepath).readlines():
    b.add(line.replace("\n", ""))
  return b

def get_plant_default_dict(species_tree):
  d = {}
  errored = False
  tree = Tree(species_tree, format = 1)
  reverse_d = {}
  for name in tree.get_leaf_names():
    sp = name.split("_")
    if (len(sp) > 2):
      new_name = "_".join(sp[0:2])
      if (new_name in reverse_d):
        print("WARNING: duplicate name between " + name + " and " + reverse_d[new_name])
        errored = True
      else:
          d[name] = "_".join(sp[0:2])
          reverse_d[new_name] = name
  #if  (errored):
  #  sys.exit(1)
  return d

def remove_species_tree_inplace(species_tree, to_remove):
  tree = Tree(species_tree, format = 1)
  to_keep = []
  for leaf in tree.get_leaves():
    if (not leaf.name in to_remove):
      to_keep.append(leaf.name)
  tree.prune(to_keep, True)
  tree.write(outfile = species_tree)
  print("UPDATE " + species_tree)

def rename_species_tree_inplace(species_tree, to_replace_1, to_replace_2):
  tree = Tree(species_tree, format = 1)
  for leaf in tree.get_leaves():
    if (leaf.name in to_replace_1):
      leaf.name = to_replace_1[leaf.name]
    elif (leaf.name in to_replace_2):
      leaf.name = to_replace_2[leaf.name]
  print("UPDATE " + species_tree)
  tree.write(outfile = species_tree)

def edit_species_tree(database, species_tree, edited_species_tree):
  tree = Tree(species_tree, format = 1)
  tree.write(outfile = edited_species_tree)
  if (database == "plants"):
    bl = load_species_blacklist(exp.ensembl_plants_species_blacklist)
    remove_species_tree_inplace(edited_species_tree, bl)
    d1 = load_species_dict(exp.ensembl_plants_species_dict)
    d2 = get_plant_default_dict(edited_species_tree)
    rename_species_tree_inplace(edited_species_tree, d1, d2)
  elif (database != "vertebrates"):
    print("Invalid database value")
    sys.exit(1)
    remove_species_tree_inplace(edited_species_tree, species_to_remove)
    rename_species_tree_inplace(edited_species_tree, species_to_replace, species_to_replace)
  


def mycmp(a, b):
  return (a > b) - (a < b) 

def gene_name(name):
  return name.replace(":", "")

class SeqEntry():
  def __init__(self, line, family):
    split = line.split(" ")
    self.species = title_node_names.get_title(split[1])
    self.gene = gene_name(split[2])
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
    return mycmp(entry1.species, entry2.species)
  if (entry1.chromozome != entry2.chromozome):
    return mycmp(entry1.chromozome, entry2.chromozome)
  return mycmp(entry1.begin, entry2.begin)

# read an ensembl tree, prune the genes that we do not consider,
# and return true if the resulting tree has more than 3 taxa
def read_ensembl_tree(line, family, trees_dict, seq_entries_dict, alignments_dico):
  tree = Tree(line.replace("\n", ""), format = 1)
  leaves = tree.get_leaves()
  filtered_leaves = []
  for leaf in leaves:
    leaf.name = gene_name(leaf.name)
    if leaf.name in seq_entries_dict and leaf.name in alignments_dico:
      filtered_leaves.append(leaf.name)
  try:
    tree.prune(filtered_leaves)
  except:
    return False
  print(len(tree.get_leaves()))
  if (len(tree.get_leaves()) > 3):
    trees_dict[family] = tree.write().replace("\n", "")
    return True
  return False


def parse_nhx_emf(emf_file, species_dict, max_families, alignments_dico):
  seq_entries_dict = {}  
  current_entries_dict = {}
  trees_dict = {}
  family_index = 0
  family = "family_" + str(family_index)
  last_was_data = False
  for line in open(emf_file).readlines():
    if (last_was_data):
      last_was_data = False
      if (read_ensembl_tree(line, family, trees_dict, current_entries_dict, alignments_dico)):
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


def parse_fasta(fasta_file):
  alignments_dico = {}
  cur_name = ""
  cur_seq = ""
  for line in open(fasta_file).readlines():
    if (line.strip() == "" or line[0] == "/"):
      continue
    if (line[0] == ">"):
      if (len(cur_name)):
        alignments_dico[cur_name] = cur_seq
      cur_name = gene_name(line[1:-1])
      cur_seq = line
    else:
      cur_seq += line
  alignments_dico[cur_name] = cur_seq
  print("number of DNA seq: " + str(len(alignments_dico)))
  return alignments_dico

def parallelized_function(params):
  datadir, family = params
  create_random_tree.create_random_tree(fam.get_alignment(datadir, family), fam.get_true_tree(datadir, family))

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


def extract_adjacencies(seq_entries_dict, datadir):
  seq_entries_list = []
  for gene in seq_entries_dict:
    seq_entries_list.append(seq_entries_dict[gene])
  seq_entries_list.sort(key = functools.cmp_to_key(compare))
  with open(fam.get_adjacencies(datadir), "w") as writer:
    for gene1, gene2 in pairwise(seq_entries_list):
      if (gene1.species == gene2.species and gene1.chromozome == gene2.chromozome):
        writer.write(gene1.gene + " " + gene2.gene + "\n")
  with open(fam.get_prefixed_adjacencies(datadir), "w") as writer:
    for gene1, gene2 in pairwise(seq_entries_list):
      if (gene1.species == gene2.species and gene1.chromozome == gene2.chromozome):
        writer.write(gene1.species + "_" + gene1.gene + " " + gene2.species + "_" + gene2.gene + "\n")

def extract_deco_mappings(seq_entries_dict, output_file):
  with open(output_file, "w") as writer:
    for gene in seq_entries_dict:
      writer.write(seq_entries_dict[gene].species + " " + gene + "\n")

def export(species_tree, seq_entries_dict, trees_dict, alignments_dico, datadir):
  fam.init_top_directories(datadir)
  per_family_seq_entries = {}
  for gene in seq_entries_dict:
    seq_entry = seq_entries_dict[gene]
    if (not seq_entry.family in per_family_seq_entries):
      per_family_seq_entries[seq_entry.family] = []
    per_family_seq_entries[seq_entry.family].append(seq_entry)
    fam.init_family_directories(datadir, seq_entry.family)

  shutil.copy(species_tree, fam.get_species_tree(datadir))
  
  families_dir = fam.get_families_dir(datadir)
  print("Number of families: " + str(len(per_family_seq_entries)))
  processed_families = 0
  for family in per_family_seq_entries:
    family_dir = fam.get_family_path(datadir, family)
    seq_entries = per_family_seq_entries[family]
    export_msa(seq_entries, alignments_dico, fam.get_alignment(datadir, family)) 
    export_mappings(seq_entries, fam.get_mappings(datadir, family))
    processed_families += 1
    with open(fam.get_true_tree(datadir, family), "w") as writer:
      writer.write(trees_dict[family])
  extract_adjacencies(seq_entries_dict, datadir)
  extract_deco_mappings(seq_entries_dict, fam.get_deco_mappings(datadir))
  fam.postprocess_datadir(datadir)

def get_species_dict(species_tree_file):
  res = {}
  tree = Tree(species_tree_file, 1)
  for leaf in tree.get_leaves():
    split = leaf.name.split("_")
    res[leaf.name] = True
  return res


def extract_from_ensembl(database, nhx_emf_file, fasta_file, species_tree, datadir, max_families):
  edited_species_tree = species_tree + ".edited"
  edit_species_tree(database, species_tree, edited_species_tree)
  new_species_tree = species_tree + ".titled"
  title_node_names.title_node_names(edited_species_tree, new_species_tree)
  species_tree = new_species_tree
  species_dict = get_species_dict(species_tree)
  print(species_dict)
  alignments_dico = parse_fasta(fasta_file)
  seq_entries_dict, trees_dict = parse_nhx_emf(nhx_emf_file, species_dict, max_families, alignments_dico)
  
  export(species_tree, seq_entries_dict, trees_dict, alignments_dico, datadir)


def get_max_families(argv):
  max_families = 10000000000
  for i in range(0, len(argv)):
    if (argv[i] == "--max-families"):
      max_families = int(argv[i+1])
      print("reading max families " + str(max_families))
  return max_families

if (__name__ == "__main__"): 
  if (len(sys.argv) < 6): 
    print("Syntax: python " + os.path.basename(__file__) + " database nhx_emf fasta species_tree output")
    print("Database can be either plants or vertebrates")
    exit(1)
  database = sys.argv[1]
  nhx_emf_file = sys.argv[2]
  fasta_file = sys.argv[3]
  species_tree = sys.argv[4]
  datadir = sys.argv[5]
  max_families = get_max_families(sys.argv)
  extract_from_ensembl(database, nhx_emf_file, fasta_file, species_tree, datadir, max_families)
  

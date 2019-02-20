import os
import sys
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "mappings"))
sys.path.insert(0, os.path.join("tools", "trees"))
sys.path.insert(0, os.path.join("tools", "msa_edition"))

import experiments as exp
import make_binary
import nhx_to_newick
import cut_node_names
import shutil
import fasta_to_mapping
import phyldog_to_treerecs_map
import align
from ete3 import Tree
from ete3 import SeqGroup
from sets import Set

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None

def get_leaf_names(tree_file):
  tree = read_tree(tree_file)
  leaves = tree.get_leaves()
  leaf_names = []
  for leaf in leaves:
    leaf_names.append(leaf.name)
  return Set(leaf_names)

def copy_remove_carriage(input_file, output_file):
  s = open(input_file).read()
  with open(output_file, "w") as writer:
    writer.write(s.replace('\r',''))

def extract_species_tree(swiss_dataset, output_dataset):
  nhx_species_tree = os.path.join(swiss_dataset, "speciestree.nhx")
  unresolved_species_tree = os.path.join(output_dataset, "unresolvedSpeciesTree.newick")
  species_tree = os.path.join(output_dataset, "speciesTree.newick")
  nhx_to_newick.nhx_to_newick(nhx_species_tree, unresolved_species_tree)
  make_binary.make_binary(unresolved_species_tree, species_tree)
  cut_node_names.cut_node_names(species_tree, species_tree, True)
  return species_tree

def remove_fasta_comments(swiss_alignment, output_alignment):
  lines = open(swiss_alignment).readlines()
  with open(output_alignment, "w") as writer:
    ok = False
    for line in lines:
      if (line.startswith(">")):
        ok = True
        sp = line.split("|")
        if (len(sp) != 1):
          line = ">" +  sp[-1].split(" ")[0]
          print(line) 
      if (ok):
        writer.write(line + "\n")

def extract_sequence(swiss_alignment, output_alignment, species_leaves):
  remove_fasta_comments(swiss_alignment, output_alignment)
  seq = SeqGroup(output_alignment)
  new_seq = SeqGroup()
  for entry in seq.get_entries():
    name = entry[0]
    species = name.split("_")[1]
    if (not species in species_leaves):
      # get rid of genes mapped to species
      # that are not in the species tree
      continue
    new_seq.set_seq(entry[0], entry[1])
  new_seq.write("fasta", output_alignment)

def extract_family(swiss_dataset, family, species_tree, output_dataset):
  print("Treating " + family)
  swiss_alignment = os.path.join(swiss_dataset, family, "sequences.fst")
  if (not os.path.isfile(swiss_alignment)):
    swiss_alignment = os.path.join(swiss_dataset, family, "sequences.txt")
  if (not os.path.isfile(swiss_alignment)):
    print("Cannot extract family " + family)
    return

  swiss_sequence_id = os.path.join(swiss_dataset, family, "sequence_identifiers.txt")
  swiss_consensus_tree = os.path.join(swiss_dataset, family, "consensus_tree.nhx")
  
  alignments_path = os.path.join(output_dataset, "alignments")
  family_path = os.path.join(output_dataset, "families", family)
  family_alignment = os.path.join(family_path, "alignment.msa")
  family_alignment_unaligned = os.path.join(family_path, "unaligned_alignment.msa")
  pargenes_alignment = os.path.join(alignments_path, family + ".fasta")
  family_mapping = os.path.join(family_path, "mapping.link")
  family_treerecsc_mapping = os.path.join(family_path, "treerecs_mapping.link")
  true_tree = os.path.join(family_path, "trueGeneTree.newick")
  
  # species tree
  exp.try_make_dir(family_path)  
  shutil.copyfile(species_tree, os.path.join(family_path, "speciesTree.newick"))
  species_leaves = get_leaf_names(species_tree)
  print(species_leaves)

  # true tree
  nhx_to_newick.nhx_to_newick(swiss_consensus_tree, true_tree)
  lines = open(true_tree).readlines()
  open(true_tree, "w").write(lines[0])
  cut_node_names.cut_node_names(true_tree, true_tree, False)

  # alignment
  extract_sequence(swiss_alignment, family_alignment_unaligned, species_leaves)
  align.align(family_alignment_unaligned, family_alignment, "MAFFT")
  shutil.copyfile(family_alignment, pargenes_alignment)
  
  # mapping
  fasta_to_mapping.fasta_to_mapping(family_alignment, family_mapping)
  phyldog_to_treerecs_map.convert(family_mapping, family_treerecsc_mapping)


def swiss_to_family(swiss_dataset, output_dataset):
  species_tree = extract_species_tree(swiss_dataset, output_dataset)
  alignments_path = os.path.join(output_dataset, "alignments")
  exp.try_make_dir(alignments_path)
  for family in os.listdir(swiss_dataset):
    if (family.startswith("ST0") and family != "ST012" and family != "ST007"):
      extract_family(swiss_dataset, family, species_tree, output_dataset)
 

swiss_dataset = os.path.join(exp.benoit_datasets_root, "families/reference_genes/original_files")
output_dataset = os.path.join(exp.benoit_datasets_root, "families/reference_genes/")

swiss_to_family(swiss_dataset, output_dataset)


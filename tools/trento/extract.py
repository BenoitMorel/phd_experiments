import shutil
import os
import sys
sys.path.insert(0, os.path.join("tools", "families"))
import fam

def extract(src, dest, families, extract_alignments = True, extract_mappings = False, extract_trees = False, gene_method = "raxml-ng", subst_model = "GTR+G"):
  os.makedirs(dest)
  tree_dir = os.path.join(dest, "trees")
  ali_dir = os.path.join(dest, "alignments")
  mappings_dir = os.path.join(dest, "mappings")
  if (extract_alignments):
    os.mkdir(ali_dir)
  if (extract_mappings):
    os.mkdir(mappings_dir)
  if (extract_trees):
    os.mkdir(tree_dir)
  index = 0
  for family in families:
    new_family = "fam_" + str(index)
    index += 1
    print(index)
    if (extract_alignments):
      alisrc = fam.get_alignment(src, family)
      alidest = os.path.join(ali_dir, new_family) + ".fasta"
      shutil.copyfile(alisrc, alidest)
    if (extract_mappings):
      mapsrc = fam.get_mappings(src, family)
      mapdest = os.path.join(mappings_dir, new_family) + ".map"
      shutil.copyfile(mapsrc, mapdest)
    if (extract_trees):
      treesrc = fam.get_gene_tree_path(src, family, gene_method, subst_model) 
      treedest = os.path.join(tree_dir, new_family) + ".newick"
      shutil.copyfile(treesrc, treedest)




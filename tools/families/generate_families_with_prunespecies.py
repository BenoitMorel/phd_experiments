import sys
import os
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
from ete3 import Tree
from ete3 import SeqGroup
import get_dico
from read_tree import read_tree
from read_tree import read_trees_list

def get_species_to_keep(inputdir, species_to_prune):
  species_tree = Tree(fam.get_species_tree(inputdir), format = 1)
  all_species = set(species_tree.get_leaf_names())
  return all_species - set(species_to_prune)



"""
  detach the current node from its parent
  return the pruned subtree (which needs to be reattached)
  return None if there is the whole subtree should be pruned
"""
def prune_rec(node, to_keep):
  node.detach()
  children = node.get_children()
  if (len(children) == 0):
    #leaf
    if (node.name in to_keep):
      return node
    else:
      return None
  ok_children = []
  for child in children:
    ok_child = prune_rec(child, to_keep)
    if (ok_child != None):
      ok_children.append(ok_child)
  if (len(ok_children) == 0):
    return None
  if (len(ok_children) == 1):
    return ok_children[0]
  
  for child in ok_children:
    node.add_child(child)

  return node


def prune(tree, to_keep):
  return prune_rec(tree, to_keep)




def extract_species_tree(inputdir, outputdir, species_to_keep):
  input_tree = fam.get_species_tree(inputdir)
  output_tree = fam.get_species_tree(outputdir)
  with open(output_tree, "w") as writer:
    for tree in read_trees_list(input_tree): 
      tree = prune(tree, species_to_keep)
      writer.write(tree.write())
      writer.write("\n")

def extract_family(family, inputdir, outputdir, gene_method, subst_model, species_to_keep):
  print(family)
  #tree = Tree(fam.get_true_tree(inputdir, family), format = 1)
  species_to_genes = get_dico.get_species_to_genes_family(inputdir, family)
  new_species_to_genes = {}
  genes_to_keep = set()

  for species in species_to_genes:
    if (species in species_to_keep):
      new_species_to_genes[species] = species_to_genes[species]
      genes_to_keep |= set(species_to_genes[species])
  if (len(genes_to_keep) < 4):
    return
  fam.init_family_directories(outputdir, family)
  fam.write_phyldog_mapping(new_species_to_genes, fam.get_mappings(outputdir, family))
  input_gene_tree = fam.build_gene_tree_path(inputdir, subst_model, family, gene_method)
  if (os.path.isfile(input_gene_tree)):
    trees = read_trees_list(input_gene_tree)
    output_gene_tree = fam.build_gene_tree_path(outputdir, subst_model, family, gene_method)
    with open(output_gene_tree, "w") as writer:
      for tree in trees:
        tree = prune(tree, genes_to_keep)
        writer.write(tree.write())
        writer.write("\n")
  input_ali = fam.get_alignment(inputdir, family)
  if (os.path.isfile(input_ali)):
    msa = SeqGroup(input_ali)
    new_msa = SeqGroup()
    for entry in msa:
      if (entry[0] in genes_to_keep):
        new_msa.set_seq(entry[0], entry[1])
    new_msa.write(outfile = fam.get_alignment(outputdir, family))

def generate(inputdir, outputdir, gene_method, subst_model, keep, species_to_prune):
  print("Species to keep/prune: " + " ".join(species_to_prune))
  fam.init_top_directories(outputdir)   
  with open(os.path.join(fam.get_misc_dir(outputdir), "info.txt"), "w") as writer:
    writer.write("Extracted from " + os.path.basename(inputdir))
    if (keep):
      writer.write(" by keeping the species:\n" + "\n".join(species_to_prune))
    else:
      writer.write(" by removing the species:\n" + "\n".join(species_to_prune))
  species_to_keep = None
  if (keep):
    species_to_keep = species_to_prune
  else:
    species_to_keep = get_species_to_keep(inputdir, species_to_prune)
  print("Species to keep: " + " ".join(species_to_keep))
  extract_species_tree(inputdir, outputdir, species_to_keep)
  families = fam.get_families_list(inputdir)
  index = 0
  for family in families:
    print("treating " + family + " " + str(index) + "/" + str(len(families)))
    extract_family(family, inputdir, outputdir, gene_method, subst_model, species_to_keep)
    index += 1
  fam.postprocess_datadir(outputdir)
  print("Result datadir in " + outputdir)
  print("Species to keep: " + " ".join(species_to_keep))


if (__name__ == "__main__"): 
  if (len(sys.argv) < 6): 
    print("Syntax: python " + os.path.basename(__file__) + " input output gene_method subst_model keep  use_species_tree species1 [species2 species3 ...]")
    print("Set keep to 0 for pruning the species, and to 1 for keeping the species")
    exit(1)
  inputdir = sys.argv[1]
  outputdir = sys.argv[2]
  gene_method = sys.argv[3]
  subst_model = sys.argv[4]
  keep = int(sys.argv[5]) != 0
  use_species_tree = int(sys.argv[6]) != 0
  species_to_prune = sys.argv[7:]
  if (use_species_tree):
    species_tree = sys.argv[7]
    species_to_prune = read_tree(species_tree).get_leaf_names()
  generate(inputdir, outputdir, gene_method, subst_model, keep, species_to_prune)




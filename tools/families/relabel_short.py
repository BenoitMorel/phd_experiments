import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "mappings"))
import fam
import ete3
import get_dico

def build_short_names_seq():
  alphabet = " ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  seq = []
  for a in alphabet:
    for b in alphabet[1:]:
      seq.append((a + b).replace(" ", ""))
  return seq

def get_short_mapping(datadir):
  true_species_tree = ete3.Tree(fam.get_species_tree(datadir), format = 1)
  short_mapping = {}
  index = 1
  short_names_seq = build_short_names_seq()
  short_mapping_file = os.path.join(fam.get_misc_dir(datadir), "short_mappings.txt")
  writer = open(short_mapping_file, "w")
  print("Writting short mapping in " + short_mapping_file)
  for name in true_species_tree.get_leaf_names():
    writer.write(short_names_seq[index] + " " + name + "\n")
    short_mapping[name] = short_names_seq[index]
    index += 1
  return short_mapping

def relabel_gene_trees(datadir, gene_trees, subst_model):
  short_mapping = get_short_mapping(datadir)
  for family in fam.get_families_list(datadir):
    try:
      gene_tree_path = fam.build_gene_tree_path(datadir, subst_model, family, gene_trees)
      output_gene_tree_path = gene_tree_path + ".shortlabels"
      tree = ete3.Tree(gene_tree_path)
      for leaf in tree.get_leaves():
        mapping = get_dico.get_gene_to_species(datadir, family)
        leaf.name = short_mapping[mapping[leaf.name]]
      open(output_gene_tree_path, "w").write(tree.write())
    except:
      print("Error for family " + family)


def relabel_species_tree(datadir, species_tree, subst_model):
  short_mapping = get_short_mapping(datadir)
  tree_path = fam.get_species_tree(datadir, subst_model, species_tree)
  print("Relabelling " + tree_path)
  tree = ete3.Tree(tree_path, format = 1)
  for leaf in tree.get_leaves():
    leaf.name = short_mapping[leaf.name]
  output_tree_path = tree_path + ".shortlabels"
  open(output_tree_path, "w").write(tree.write())
  print(tree)
  print("Writting new tree in " + output_tree_path)

if (__name__== "__main__"):
  if (len(sys.argv) == 5):
    mode = sys.argv[1]
    datadir = sys.argv[2]
    method = sys.argv[3]
    subst_model = sys.argv[4]
    if (mode == "genes"):
      relabel_gene_trees(datadir, method, subst_model)
    elif(mode == "species"):
      relabel_species_tree(datadir, method, subst_model)
    else:
      print("Invalid mode " + mode)
  else:
    print("Syntax:")
    print(" mode datadir tree_method subst_model")
    sys.exit(0)








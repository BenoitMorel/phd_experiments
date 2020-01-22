import sys
import os
import shutil
import random
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
import ete3
import get_dico

def cp(src, dest):
  shutil.copy(src, dest)

def get_gene_list(input_datadir, family):
  res = []
  for mapping in get_dico.get_gene_to_species(input_datadir, family):
    res.append(mapping)
  return res

def get_sampled_gene_set(input_gene_list, random_sampling_ratio):
  res = []
  for gene in input_gene_list:
    if (random_sampling_ratio > random.uniform(0.0, 1.0)):
      res.append(gene)
  return set(res)

def copy_mappings(input_datadir, output_datadir, family, output_gene_set):
  mappings = get_dico.get_species_to_genes_family(input_datadir, family)
  for species in mappings:
    sampled_genes_set = output_gene_set.intersection(mappings[species])
    mappings[species] = list(sampled_genes_set)
  get_dico.export_species_to_genes(mappings, output_datadir, family)

def copy_gene_tree(input_datadir, output_datadir, family, output_gene_set):
  old_gene_tree = fam.get_true_tree(input_datadir, family)
  new_gene_tree = fam.get_true_tree(output_datadir, family)
  tree = ete3.Tree(old_gene_tree, format=1)
  tree.prune(output_gene_set, preserve_branch_length = True)
  open(new_gene_tree, "w").write(tree.write())

def copy_alignment(input_datadir, output_datadir, family, output_gene_set):
  old_alignment = ete3.SeqGroup(fam.get_alignment(input_datadir, family), "fasta")
  new_alignment = ete3.SeqGroup()
  for entry in old_alignment.iter_entries():
    if (entry[0] in output_gene_set):
      new_alignment.set_seq(entry[0], entry[1])
  new_alignment.write("fasta", fam.get_alignment(output_datadir, family))

def copy_sampled_genes(input_datadir, output_datadir, family, output_gene_set):
  fam.init_family_directories(output_datadir, family)
  copy_mappings(input_datadir, output_datadir, family, output_gene_set)
  copy_gene_tree(input_datadir, output_datadir, family, output_gene_set)
  copy_alignment(input_datadir, output_datadir, family, output_gene_set)


def sample_family(input_datadir, output_datadir, family, random_sampling_ratio):
  input_gene_list = get_gene_list(input_datadir, family)
  output_gene_set = get_sampled_gene_set(input_gene_list, random_sampling_ratio)
  print(str(len(input_gene_list)) + " " + str(len(output_gene_set)))
  if (len(output_gene_set) > 3):
    copy_sampled_genes(input_datadir, output_datadir, family, output_gene_set)

def sample_missing_data(input_datadir, output_datadir, random_sampling_ratio):
  random.seed = 41
  fam.init_top_directories(output_datadir)
  families = fam.get_families_list(input_datadir)
  cp(fam.get_species_tree(input_datadir), fam.get_species_tree(output_datadir))
  for family in families:
    sample_family(input_datadir, output_datadir, family, random_sampling_ratio)
  fam.postprocess_datadir(output_datadir)

if (__name__ == "__main__"): 
  if (len(sys.argv) < 4): 
    print("Syntax: python " + os.path.basename(__file__) + " input_datadir output_datadir random_sampling_ratio")
  input_datadir = sys.argv[1]
  output_datadir = sys.argv[2]
  random_sampling_ratio = float(sys.argv[3])
  sample_missing_data(input_datadir, output_datadir, random_sampling_ratio)



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

def get_sampled_gene_set(input_gene_list, gene_to_species, species_missing_ratio):
  res = []
  for gene in input_gene_list:
    missing_ratio = species_missing_ratio[gene_to_species[gene]]
    if (missing_ratio > random.uniform(0.0, 1.0)):
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


def sample_family(input_datadir, output_datadir, family, species_missing_ratio):
  input_gene_list = get_gene_list(input_datadir, family)
  gene_to_species = get_dico.get_gene_to_species(input_datadir, family)
  output_gene_set = get_sampled_gene_set(input_gene_list, gene_to_species, species_missing_ratio)
  if (len(output_gene_set) > 3):
    copy_sampled_genes(input_datadir, output_datadir, family, output_gene_set)

def get_species_missing_ratio(input_datadir, mu, theta):
  species_tree = ete3.Tree(fam.get_species_tree(input_datadir))
  species_missing_ratio = {}
  for species in species_tree.get_leaf_names():
    v = random.normalvariate(mu, theta)
    species_missing_ratio[species] = max(0, min(1.0, v))
  return species_missing_ratio
  
def export_missing_ratios(output_datadir, species_missing_ratio):
  print("Missing data ratios:")
  with open(fam.get_missing_data_file(output_datadir), "w") as writer:
    for k,v in sorted(species_missing_ratio.items(), key=lambda x: x[1]):
      writer.write(k + " " + str(v) + "\n")
      print(k + " " + str(v))
      
def sample_missing_data(input_datadir, output_datadir, mu, theta):
  random.seed = 41
  if (os.path.exists(output_datadir)):
    print("PATH " + output_datadir + " already exists!!!")
    sys.exit(1)
  fam.init_top_directories(output_datadir)
  families = fam.get_families_list(input_datadir)
  cp(fam.get_species_tree(input_datadir), fam.get_species_tree(output_datadir))
  cp(fam.get_discordance_file(input_datadir), fam.get_discordance_file(output_datadir))
  species_missing_ratio = get_species_missing_ratio(input_datadir, mu, theta)
  for family in families:
    sample_family(input_datadir, output_datadir, family, species_missing_ratio)
  export_missing_ratios(output_datadir, species_missing_ratio)
  fam.postprocess_datadir(output_datadir)

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5): 
    print("Syntax: python " + os.path.basename(__file__) + " input_datadir output_datadir mu theta")
  input_datadir = sys.argv[1]
  output_datadir = sys.argv[2]
  mu = float(sys.argv[3])
  theta = float(sys.argv[4])
  sample_missing_data(input_datadir, output_datadir, mu, theta)



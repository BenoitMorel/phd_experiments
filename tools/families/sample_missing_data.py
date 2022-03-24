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

def get_sampled_gene_set(input_gene_list, gene_to_species, species_sample_prob, fam_miss_prob):
  res = []
  for gene in input_gene_list:
    miss_species = species_sample_prob[gene_to_species[gene]]
    sample_prob = (1.0 - miss_species) * (1.0 - fam_miss_prob)
    if (sample_prob > random.uniform(0.0, 1.0)):
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

def beta_reparametrized(mean, shape):
  if (mean == 1.0):
    return 1.0
  if (mean == 0.0):
    return 0.0
  alpha = mean * shape
  beta = (1 - mean) * shape
  return random.betavariate(alpha, beta)



def sample_family_fam(input_datadir, output_datadir, family, species_sample_prob, miss_fam, sampled_species):
  input_gene_list = get_gene_list(input_datadir, family)
  gene_to_species = get_dico.get_gene_to_species(input_datadir, family)
  fam_miss_prob = beta_reparametrized(miss_fam, 5.0) 
  output_gene_set = get_sampled_gene_set(input_gene_list, gene_to_species, species_sample_prob, fam_miss_prob)
  if (len(output_gene_set) > 3):
    copy_sampled_genes(input_datadir, output_datadir, family, output_gene_set)
    for gene in output_gene_set:
      sampled_species.add(gene_to_species[gene])


def get_species_sample_prob(input_datadir, mu, theta = 5.0):
  species_tree = ete3.Tree(fam.get_species_tree(input_datadir))
  species_sample_prob = {}
  for species in species_tree.get_leaf_names():
    v = beta_reparametrized(mu, theta)
    species_sample_prob[species] = v # max(0.05, v)
  return species_sample_prob
  
def export_sample_probs(output_datadir, species_sample_prob):
  print("Missing data ratios:")
  with open(fam.get_missing_data_file(output_datadir), "w") as writer:
    for k,v in sorted(species_sample_prob.items(), key=lambda x: x[1]):
      writer.write(k + " " + str(v) + "\n")
      
def sample_missing_data(input_datadir, output_datadir, miss_species, miss_fam):
  random.seed = 42
  if (os.path.exists(output_datadir)):
    print("PATH " + output_datadir + " already exists!!!")
    sys.exit(1)
  fam.init_top_directories(output_datadir)
  families = fam.get_families_list(input_datadir)
  cp(fam.get_discordance_file(input_datadir), fam.get_discordance_file(output_datadir))
  species_sample_prob = get_species_sample_prob(input_datadir, miss_species)
  sampled_species = set()
  for family in families:
    sample_family_fam(input_datadir, output_datadir, family, species_sample_prob, miss_fam, sampled_species)
  export_sample_probs(output_datadir, species_sample_prob)
  fam.postprocess_datadir(output_datadir)
  species_tree = ete3.Tree(fam.get_species_tree(input_datadir), format = 1)
  species_tree.prune(sampled_species)
  species_tree.write(format = 1, outfile = fam.get_species_tree(output_datadir), dist_formatter="%0.8f")

if (__name__ == "__main__"): 
  if (len(sys.argv) < 5): 
    print("Syntax: python " + os.path.basename(__file__) + " input_datadir output_datadir miss_species miss_fam")
  input_datadir = sys.argv[1]
  output_datadir = sys.argv[2]
  miss_species = float(sys.argv[3])
  miss_fam = float(sys.argv[4])
  sample_missing_data(input_datadir, output_datadir, miss_species, miss_fam)



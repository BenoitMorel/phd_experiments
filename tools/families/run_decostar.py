import sys
import os
import subprocess
import shutil
import saved_metrics
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/notung')
import experiments as exp
import fam
from ete3 import Tree


"""
  load the gene to species dictionary, to prefix gene names (see get_prefixed_name)
"""
def load_mapping(datadir):
  mapping = {}
  for line in open(fam.get_deco_mappings(datadir)).readlines():
    s = line.split()
    mapping[s[1]] = s[0]
  return mapping

"""
  decostar: gene names must be prefixed with the species name
"""
def get_prefixed_name(gene_name, mapping):
  return mapping[gene_name] + "_" + gene_name

def prefix_tree(tree_file, prefixed_tree_file, mapping):
  tree = Tree(open(tree_file).read(), format=1)
  for leaf in tree.get_leaves():
    leaf.name = get_prefixed_name(leaf.name, mapping)
  open(prefixed_tree_file, "w").write(tree.write())

def get_decostar_run_dir(datadir, method):
  return os.path.join(datadir, "runs", "decostar_run", method)

"""
  decostar takes as input a file with all the gene trees (with prefixed gene names)
  build this file and return a path to it
"""
def build_trees_file(datadir, method, mapping):
  run_dir = get_decostar_run_dir(datadir, method)
  prefixed_trees_dir = os.path.join(run_dir, "prefixed_trees")
  os.mkdir(prefixed_trees_dir)
  trees_file = os.path.join(get_decostar_run_dir(datadir, method), "tree_files.txt")
  with open(trees_file, "w") as writer:
    for family in fam.getFamiliesList(datadir):
      tree_file = fam.get_gene_tree(fam.getFamily(datadir, family), method)
      prefixed_tree_file = os.path.join(prefixed_trees_dir, family + ".prefixed.newick")
      prefix_tree(tree_file, prefixed_tree_file, mapping)
      writer.write(prefixed_tree_file)
      writer.write("\n")
  return trees_file

def build_decostar_conf(datadir, method, trees_file, mapping):
  run_dir = get_decostar_run_dir(datadir, method)
  conf_file = os.path.join(run_dir, "conf.txt")
  with open(conf_file, "w") as writer:
    writer.write("gene.distribution.file=" + trees_file + "\n")
    #writer.write("genes_file " + fam.get_deco_mappings(datadir) + "\n")
    writer.write("species.file=" + fam.get_species_tree(datadir) + "\n")
    writer.write("adjacencies.file=" + fam.get_prefixed_adjacencies(datadir) + "\n")
    writer.write("dated.species.tree=0" + "\n")
    writer.write("with.transfer=0" + "\n")
    writer.write("already.reconciled=0" + "\n")
    writer.write("output.prefix=decostar_run" + "\n")
    writer.write("output.dir=" + run_dir + "\n")
  return os.path.abspath(conf_file)


def run_decostar(datadir, method, force):
  decostar_run_dir =  get_decostar_run_dir(datadir, method)
  if (not force and os.path.isdir(decostar_run_dir)):
    print("directory already exists")
    return
  shutil.rmtree(decostar_run_dir, True)
  os.makedirs(decostar_run_dir)
  mapping = load_mapping(datadir)
  logs_file = os.path.join(decostar_run_dir, "logs.txt")
  trees_file = build_trees_file(datadir, method, mapping)
  conf_file = build_decostar_conf(datadir, method, trees_file, mapping)
  command = [exp.decostar_exec, "parameter.file=" + conf_file]
  print("Executing " + " ".join(command))
  print("Output in " + logs_file)
  subprocess.check_call(command, stdout = open(logs_file, "w"))
  print("Results in : " + decostar_run_dir)

def get_distance_to_1(count_vector):
  distance = 0.0
  for i in range(2, len(count_vector)):
    distance += abs(float(i) - 1.0) * float(count_vector[i])
  return distance 

def get_distance(per_species_adjacencies_left, per_species_adjacencies_right, species_to_analyze, all_genes, norm):
  pass


"""
  Analyze decostar outputs: compute, print and save all "ancestral metrics" 
  (adjacencies and gene counts)
"""
def analyze_decostar_output(datadir, method):
  decostar_run_dir =  get_decostar_run_dir(datadir, method)
  adjacencies_file = os.path.join(decostar_run_dir, "decostar_run.adjacencies.txt")
  # per_species_adjacencies_left[s][g] = number left neighbors of gene g in species s
  per_species_adjacencies_left = {}
  # per_species_adjacencies_right[s][g] = number right neighbors of gene g in species s
  per_species_adjacencies_right = {}
  # dictionary of all ancestral genes
  all_genes = {}
  # fill per_species_adjacencies_left, per_species_adjacencies_right and all_genes
  for line in open(adjacencies_file).readlines():
    if (not "|" in line):
      continue
    split = line.split()
    species = split[0]
    gene1 = split[1]
    gene2 = split[2]
    all_genes[gene1] = 1
    all_genes[gene2] = 1
    if (not species in per_species_adjacencies_left):
      per_species_adjacencies_left[species] = {}
      per_species_adjacencies_right[species] = {}
    species_adjacencies_left = per_species_adjacencies_left[species]
    species_adjacencies_right = per_species_adjacencies_right[species]
    if (not gene1 in species_adjacencies_left):
      species_adjacencies_right[gene1] = 0
      species_adjacencies_left[gene1] = 0
    species_adjacencies_left[gene1] += 1
    if (not gene2 in species_adjacencies_right):
      species_adjacencies_left[gene2] = 0
      species_adjacencies_right[gene2] = 0
    species_adjacencies_right[gene2] += 1
  
  # compute sorted species vector 
  sorted_species_list = []
  for species in per_species_adjacencies_left:
    sorted_species_list.append(int(species))
  sorted_species_list.sort()
 
  # root only
  #species_to_analyze = sorted_species_list[-1:]
  # all branches
  species_to_analyze = sorted_species_list

  
  # get max adj
  max_adj_count = 0
  total_adj_count = 0
  for species in species_to_analyze:
    for gene in per_species_adjacencies_left[str(species)]:
      max_adj_count = max(max_adj_count, per_species_adjacencies_left[str(species)][gene])
      total_adj_count += per_species_adjacencies_left[str(species)][gene]
    for gene in per_species_adjacencies_right[str(species)]:
      max_adj_count = max(max_adj_count, per_species_adjacencies_right[str(species)][gene])
      total_adj_count += per_species_adjacencies_right[str(species)][gene]
  ancestral_genes_count = len(all_genes)
  # extract counts
  adj_count_left = [0] * (max_adj_count + 1)
  adj_count_right = [0] * (max_adj_count + 1)
  total_classes = 0
  for species in species_to_analyze:
    species = str(species)
    species_adjacencies_left = per_species_adjacencies_left[species]
    for gene in species_adjacencies_left:
      adj_count_left[species_adjacencies_left[gene]] += 1
      total_classes += 1
    species_adjacencies_right = per_species_adjacencies_right[species]
    for gene in species_adjacencies_right:
      adj_count_right[species_adjacencies_right[gene]] += 1
      total_classes += 1
  distance = get_distance_to_1(adj_count_left) + get_distance_to_1(adj_count_right)
  distance = float(distance) / float(total_classes)
  print("Total number of adjacencies: " + str(total_adj_count))
  print("Total number of adjacencies class: " + str(total_classes))
  print("Total number of genes: " + str(ancestral_genes_count))
  print("distance to 1: " + str(distance))
  saved_metrics.save_metrics(datadir, method, distance, "syntheny_score_1")
  saved_metrics.save_metrics(datadir, method, ancestral_genes_count, "gene_count")
  saved_metrics.save_metrics(datadir, method, total_adj_count, "syntheny_count")

def run_on_analized_methods(datadir):
  for method in saved_metrics.get_metrics_methods(datadir, "average_rrf"):
    run_decostar(datadir, method, force = True)
    analyze_decostar_output(datadir, method)
  run_decostar(datadir, "true", force = True)
  analyze_decostar_output(datadir, "true")

if (__name__ == "__main__"): 
  if (len(sys.argv) < 3): 
     print("Syntax: python " + os.path.basename(__file__) + " datadir method [run_decostar]")
     print("method can also be set to all")
     exit(1)
  datadir = sys.argv[1]
  method = sys.argv[2]
  force = False

  if (method == "all"):
    run_on_analized_methods(datadir)
  else:
    if (len(sys.argv) > 3):
      force = int(sys.argv[3])
    run_decostar(datadir, method, force)
    analyze_decostar_output(datadir, method)

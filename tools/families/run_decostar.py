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
  return prefixed gene names with species names, as required by decostar
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
  build the file containing the list of input trees for decostar
"""
def build_trees_file(datadir, method, mapping):
  run_dir = get_decostar_run_dir(datadir, method)
  prefixed_trees_dir = os.path.join(run_dir, "prefixed_trees")
  os.mkdir(prefixed_trees_dir)
  trees_file = os.path.join(get_decostar_run_dir(datadir, method), "tree_files.txt")
  with open(trees_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      tree_file = fam.get_gene_tree(fam.get_family_path(datadir, family), method)
      prefixed_tree_file = os.path.join(prefixed_trees_dir, family + ".prefixed.newick")
      prefix_tree(tree_file, prefixed_tree_file, mapping)
      writer.write(prefixed_tree_file)
      writer.write("\n")
  return trees_file

"""
  build decostar conf file 
"""
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

"""
  generate decostar inputs and run it
"""
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

"""
  get the sum of distances to one from a count vector
"""
def get_distance_to_1(count_vector):
  distance = 0.0
  for i in range(0, len(count_vector)):
    distance += abs(float(i) - 1.0) * float(count_vector[i])
  return distance 

def get_distance(adjacencies_left, adjacencies_right, per_species_genes, species_to_analyze, check_orientation = True, norm = 1.0):
  distance = 0.0
  for species in species_to_analyze:
    for gene in per_species_genes[species]:
      if (check_orientation):
        distance += float(abs(1.0 - float(adjacencies_left[species][gene])))
        distance += float(abs(1.0 - float(adjacencies_right[species][gene])))
      else:
        adj_count = adjacencies_left[species][gene] + adjacencies_right[species][gene]
        distance += float(abs(2.0 - float(adj_count)))
  return distance / norm

"""
  read decostart output and store adjacencies and genes into directionaries
"""
def extract_decostar_adjacencies(adjacencies_file, adjacencies_left, adjacencies_right, per_species_genes): 
  for line in open(adjacencies_file).readlines():
    if (not "|" in line):
      continue
    split = line.split()
    species = split[0]
    gene_left = split[1]
    gene_right = split[2]
    if (not species in per_species_genes):
      per_species_genes[species] = {}
      adjacencies_left[species] = {}
      adjacencies_right[species] = {}
    if (not gene_left in per_species_genes[species]):
      per_species_genes[species][gene_left] = True
      adjacencies_right[species][gene_left] = 0
      adjacencies_left[species][gene_left] = 0
    adjacencies_left[species][gene_left] += 1
    if (not gene_right in per_species_genes[species]):
      per_species_genes[species][gene_right] = True
      adjacencies_left[species][gene_right] = 0
      adjacencies_right[species][gene_right] = 0
    adjacencies_right[species][gene_right] += 1


def compute_and_save_scores(datadir, method, adjacencies_left, adjacencies_right, per_species_genes):
  # compute sorted species vector 
  sorted_species_list = []
  for species in per_species_genes:
    sorted_species_list.append(int(species))
  sorted_species_list.sort()
  sorted_species_list = map(str, sorted_species_list)

  total_adj_count = 0
  for species in sorted_species_list:
    for gene in adjacencies_left[str(species)]:
      total_adj_count += adjacencies_left[str(species)][gene]
    for gene in adjacencies_right[str(species)]:
      total_adj_count += adjacencies_right[str(species)][gene]
  ancestral_genes_count = 0
  for species in per_species_genes:
    ancestral_genes_count += len(per_species_genes[species])
  
  
  dist_oriented = get_distance(adjacencies_left, adjacencies_right, per_species_genes, sorted_species_list, True, 1.0)
  dist_unoriented = get_distance(adjacencies_left, adjacencies_right, per_species_genes, sorted_species_list, False, 1.0)
  dist_unoriented_weightgenes = get_distance(adjacencies_left, adjacencies_right, per_species_genes, sorted_species_list, False, ancestral_genes_count)
  dist_unoriented_weightadj = get_distance(adjacencies_left, adjacencies_right, per_species_genes, sorted_species_list, False, total_adj_count)

  saved_metrics.save_metrics(datadir, method, ancestral_genes_count, "total_gene_count")
  saved_metrics.save_metrics(datadir, method, total_adj_count, "total_adj_count")
  saved_metrics.save_metrics(datadir, method, dist_oriented, "adj_oriented")
  saved_metrics.save_metrics(datadir, method, dist_unoriented, "adj_unoriented")
  saved_metrics.save_metrics(datadir, method, dist_unoriented_weightgenes, "adj_unoriented_weightgenes")
  saved_metrics.save_metrics(datadir, method, dist_unoriented_weightadj, "adj_unoriented_weightadj")
  
  print("Non weighted oriented distance: " + str(dist_oriented))
  print("Non weighted unoriented distance: " + str(dist_unoriented))
  print("Gene-weighted unoriented distance: " + str(dist_unoriented_weightgenes))
  print("Adj-weighted unoriented distance: " + str(dist_unoriented_weightadj))
  print("Total number of ancestral adjacencies: " + str(total_adj_count))
  print("Total number of ancestral genes: " + str(ancestral_genes_count))
   

"""
  Analyze decostar outputs: compute, print and save all "ancestral metrics" 
  (adjacencies and gene counts)
"""
def analyze_decostar_output(datadir, method):
  print("Run decostar " + datadir + " " + method)
  decostar_run_dir =  get_decostar_run_dir(datadir, method)
  adjacencies_file = os.path.join(decostar_run_dir, "decostar_run.adjacencies.txt")
  # adjacencies_left[s][g] = number left neighbors of gene g in species s
  adjacencies_left = {}
  # adjacencies_right[s][g] = number right neighbors of gene g in species s
  adjacencies_right = {}
  # dictionary of all ancestral genes
  per_species_genes = {}
  extract_decostar_adjacencies(adjacencies_file, adjacencies_left, adjacencies_right, per_species_genes) 
  compute_and_save_scores(datadir, method, adjacencies_left, adjacencies_right, per_species_genes) 

def run_on_analized_methods(datadir):
  for method in fam.get_ran_methods(datadir):
    run_decostar(datadir, method, force = False)
    analyze_decostar_output(datadir, method)
  run_decostar(datadir, "true", force = False)
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

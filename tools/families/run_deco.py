import sys
import os
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import experiments as exp
import fam




def get_deco_run_dir(datadir, method):
  return os.path.join(datadir, "runs", "deco_run", method)

def build_trees_file(datadir, method):
  trees_file = os.path.join(get_deco_run_dir(datadir, method), "trees_file.newick")
  with open(trees_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      tree_file = fam.get_gene_tree(fam.get_family_path(datadir, family), method)
      writer.write(open(tree_file).read())
      writer.write("\n")
  return trees_file

def build_deco_conf(datadir, method, trees_file):
  run_dir = get_deco_run_dir(datadir, method)
  conf_file = os.path.join(run_dir, "conf.txt")
  with open(conf_file, "w") as writer:
    writer.write("trees_file " + trees_file + "\n")
    writer.write("genes_file " + fam.get_deco_mappings(datadir) + "\n")
    writer.write("species_file " + fam.get_species_tree(datadir) + "\n")
    writer.write("adjacencies_file " + fam.get_adjacencies(datadir) + "\n")
    writer.write("exp_name dec_run" + "\n")
    writer.write("directory " + run_dir + "\n")
    writer.write("ReconcilDone false" + "\n")
    writer.write("INPUT_FORMAT 0" + "\n")
    writer.write("OUTPUT_FORMAT 1" + "\n")
    writer.write("sep |" + "\n")
    writer.write("Adj_percentage 0" + "\n")
    writer.write("Gain 2" + "\n")
    writer.write("Break 1" + "\n")
  return conf_file


def run_deco(datadir, method):
  deco_run_dir =  get_deco_run_dir(datadir, method)
  shutil.rmtree(deco_run_dir, True)
  os.makedirs(deco_run_dir)
  logs_file = os.path.join(deco_run_dir, "logs.txt")
  trees_file = build_trees_file(datadir, method)
  conf_file = build_deco_conf(datadir, method, trees_file)
  command = [exp.deco_exec, conf_file]
  print("Executing " + " ".join(command))
  print("Output in " + logs_file)
  subprocess.check_call(command, stdout = open(logs_file, "w"))
  print("Results in : " + deco_run_dir)

def analyze_deco_output(datadir, method):
  deco_run_dir =  get_deco_run_dir(datadir, method)
  adjacencies_file = os.path.join(deco_run_dir, "dec_run_OUTPUT_adjacencies")
  per_species_adjacencies = {}
  # fill adjacencies li
  for line in open(adjacencies_file).readlines():
    split = line.split()
    species = split[0]
    gene1 = split[1]
    gene2 = split[2]
    if (not species in per_species_adjacencies):
      per_species_adjacencies[species] = {}
    species_adjacencies = per_species_adjacencies[species]
    if (not gene1 in species_adjacencies):
      species_adjacencies[gene1] = 0
    species_adjacencies[gene1] += 1
    if (not gene2 in species_adjacencies):
      species_adjacencies[gene2] = 0
    species_adjacencies[gene2] += 1
  
  # sort by species
  sorted_species_list = []
  for species in per_species_adjacencies:
    sorted_species_list.append(int(species))
  sorted_species_list.sort()
  #species_to_analyze = sorted_species_list[-1:]
  species_to_analyze = sorted_species_list
  
  # get max adj
  max_adj_count = 0
  total_adj_count = 0
  for species in species_to_analyze:
    for gene in per_species_adjacencies[str(species)]:
      max_adj_count = max(max_adj_count, per_species_adjacencies[str(species)][gene])
      total_adj_count += per_species_adjacencies[str(species)][gene]
  print("Max number of adjacencies: " + str(max_adj_count))
  print("Total number of adjacencies: " + str(total_adj_count))
  

  # extract counts
  adj_count = [0] * (max_adj_count + 1)
  total_classes = 0
  for species in species_to_analyze:
    species = str(species)
    species_adjacencies = per_species_adjacencies[species]
    for gene in species_adjacencies:
      adj_count[species_adjacencies[gene]] += 1
      total_classes += 1
  for count in adj_count[1:]:
    sys.stdout.write(" " + str(count))
  print("")
  for count in adj_count[1:]:
    sys.stdout.write(" " + "%.2f" % (float(count) / float(total_classes)))
  print("")


if (__name__ == "__main__"): 
  if (len(sys.argv) < 3): 
     print("Syntax: python " + os.path.basename(__file__) + " datadir method [run_deco]")
     exit(1)
  datadir = sys.argv[1]
  method = sys.argv[2]
  do_run_deco = True
  if (len(sys.argv) > 3):
    do_run_deco = int(sys.argv[3])
  if (do_run_deco):
    run_deco(datadir, method)
  analyze_deco_output(datadir, method)

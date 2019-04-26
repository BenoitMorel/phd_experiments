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
    for family in fam.getFamiliesList(datadir):
      tree_file = fam.get_gene_tree(fam.getFamily(datadir, family), method)
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
  trees_file = build_trees_file(datadir, method)
  conf_file = build_deco_conf(datadir, method, trees_file)
  command = [exp.deco_exec, conf_file]
  print("Executing " + " ".join(command))
  subprocess.check_call(command)
  print("Results in : " + deco_run_dir)

def analyze_deco_output(datadir, method):
  deco_run_dir =  get_deco_run_dir(datadir, method)
  adjacencies_file = os.path.join(deco_run_dir, "dec_run_OUTPUT_adjacencies")
  per_species_adjacencies = {}
  for line in open(adjacencies_file).readlines():
    split = line.split()
    species = split[0]
    gene1 = split[1]
    gene2 = split[2]
    if (not species in per_species_adjacencies):
      per_species_adjacencies[species] = {}
    species_adjacencies = per_species_adjacencies[species]
    if (not gene1 in species_adjacencies):
      species_adjacencies[gene1] = []
    species_adjacencies[gene1].append(gene2)
    if (not gene2 in species_adjacencies):
      species_adjacencies[gene2] = []
    species_adjacencies[gene2].append(gene1)
  for species in per_species_adjacencies:
    adj1 = 0
    adj2 = 0
    adj3 = 0
    species_adjacencies = per_species_adjacencies[species]
    for gene in species_adjacencies:
      if (len(species_adjacencies[gene]) == 1):
        adj1 += 1
      elif(len(species_adjacencies[gene]) == 2):
        adj2 += 1
      else:
        adj3 += 1
    print(species + " " + str(adj1) + " " + str(adj2) + " " + str(adj3))

if (__name__ == "__main__"): 
  if (len(sys.argv) < 3): 
     print("Syntax: python " + os.path.basename(__file__) + " datadir method [run_deco]")
     exit(1)
  datadir = sys.argv[1]
  method = sys.argv[2]
  run_deco = True
  if (sys.argv > 3):
    run_deco = int(sys.argv[3])
  if (run_deco):
    run_deco(datadir, method)
  analyze_deco_output(datadir, method)

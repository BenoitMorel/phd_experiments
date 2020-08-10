import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import time
import saved_metrics
import ete3
import get_dico
import random
import species_analyze
import time

def build_supermatrix(datadir, subst_model, supermatrix_path, partition_path, concatenation_mode):
  
  all_species = ete3.Tree(fam.get_species_tree(datadir), 1).get_leaf_names()
  partition_writer = open(partition_path, "w")
  offset = 1
  single = (concatenation_mode == "single")
  use_all_genes = (concatenation_mode == "max")
  treated = 0
  to_treat = len(fam.get_families_list(datadir))
  columns = {}
  for species in all_species:
    columns[species] = []
  for family in fam.get_families_list(datadir):
    if (treated % 1000 == 0):
      print("Building supermatrix: " + str(treated) + "/" + str(to_treat))
    treated += 1
    seqgroup = ete3.SeqGroup(fam.get_alignment(datadir, family))
    seq_len = len(seqgroup.get_entries()[0][1])
    gaps = "-" * seq_len
    species_to_genes = get_dico.get_species_to_genes_family(datadir, family)
    species_to_sample = {}
    skip_family = False
    for species in all_species:
      if (species in species_to_genes):
        if (single and len(species_to_genes[species]) > 1):
          skip_family = True
        random.shuffle(species_to_genes[species])
    if (skip_family):
      continue
    while (len(species_to_genes) > 3):
      for species in all_species:
        seq = gaps
        if (species in species_to_genes):
          seq = seqgroup.get_seq(species_to_genes[species][-1])
          species_to_genes[species].pop()
          if (len(species_to_genes[species]) == 0):
            del species_to_genes[species]
        columns[species].append(seq)
      partition_writer.write(subst_model + ", " + family + " = ")
      partition_writer.write(str(offset) + "-" + str(offset + seq_len - 1))
      partition_writer.write("\n")
      offset += seq_len
      if (not use_all_genes):
        break
  partition_writer.close()
  print("Writing the supermatrix...")
  supermatrix = ete3.SeqGroup()
  print("Number of partitions:" + str(len(columns[all_species[0]])))
  print("Number of sites: " + str(offset))
  for species in all_species:
    print("Add " + species + " to supermatrix...")
    supermatrix.set_seq(species, "".join(columns[species]))
  supermatrix.write("phylip_relaxed", supermatrix_path)
  print("End of writing the supermatrix")
  return offset

def run_raxml(subst_model, cores, run_dir, supermatrix_path, partition_path):
  command = []
  command.append("mpiexec")
  command.append("-np")
  command.append(str(cores))
  command.append(exp.raxml_exec)
  command.append("--search1")
  command.append("--msa")
  command.append(supermatrix_path)
  command.append("--model")
  #command.append(subst_model)
  command.append(partition_path)
  command.append("--prefix")
  command.append(os.path.join(run_dir, "concatenation"))
  command.append("--seed")
  command.append("40")
  command.append("--force")
  FNULL = open(os.devnull, 'w')
  print("running " + " ".join(command))
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
  stdout, stderr = process.communicate()

def run_concatenation(datadir, concatenation_mode,  subst_model, cores):
  run_name = "concatenation-" + concatenation_mode
  run_dir = fam.get_run_dir(datadir, subst_model,  run_name)
  shutil.rmtree(run_dir, True)
  os.makedirs(run_dir)
  supermatrix_path = os.path.join(run_dir, "supermatrix.fasta")
  partition_path = os.path.join(run_dir, "supermatrix.part")
  sites = build_supermatrix(datadir, subst_model, supermatrix_path, partition_path, concatenation_mode)
  cores = min(cores, int(sites / 500))
  start = time.time()
  
  run_raxml(subst_model, cores, run_dir, supermatrix_path, partition_path)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  raxml_tree = os.path.join(run_dir, "concatenation.raxml.bestTree")
  dest = fam.get_species_tree(datadir, subst_model, run_name)
  shutil.copy(raxml_tree, dest)

if __name__ == "__main__":
  if (len(sys.argv) != 5):
    print("syntax: python run_concatenation.py datadir concatenation_mode subst_model cores")
    print("concatenation modes can be: ")
    print("- min: randomly take ONE gene from each family and each species")
    print("- max: apply min, remove the selected genes, and restart until there is not gene left")
    print("- single: only take single-copy gene families")
    sys.exit(1)
  datadir = sys.argv[1]
  concatenation_mode = sys.argv[2]
  subst_model = sys.argv[3]
  cores = int(sys.argv[4])
  assert(concatenation_mode in ["min", "max", "single"])
  run_concatenation(datadir, concatenation_mode, subst_model, cores)
  species_analyze.analyze(datadir)

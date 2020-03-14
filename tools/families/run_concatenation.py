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

def build_supermatrix(datadir, subst_model, supermatrix_path, partition_path):
  all_species = ete3.Tree(fam.get_species_tree(datadir), 1).get_leaf_names()
  supermatrix = ete3.SeqGroup()
  partition_writer = open(partition_path, "w")
  offset = 1
  for species in all_species:
    supermatrix.set_seq(species, "")
  for family in fam.get_families_list(datadir):
    seqgroup = ete3.SeqGroup(fam.get_alignment(datadir, family))
    seq_len = len(seqgroup.get_entries()[0][1])
    gaps = "-" * seq_len
    species_to_genes = get_dico.get_species_to_genes_family(datadir, family)
    species_to_sample = {}
    for species in all_species:
      seq = gaps
      if (species in species_to_genes):
        seq = seqgroup.get_seq(random.choice(species_to_genes[species]))
      supermatrix.set_seq(species, supermatrix.get_seq(species) + seq)
    partition_writer.write(subst_model + ", " + family + " = ")
    partition_writer.write(str(offset) + "-" + str(offset + seq_len - 1))
    partition_writer.write("\n")
    offset += seq_len
  supermatrix.write("fasta", supermatrix_path)
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
  command.append("43")
  command.append("--force")
  FNULL = open(os.devnull, 'w')
  print("running " + " ".join(command))
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
  stdout, stderr = process.communicate()

def run_concatenation(datadir, subst_model, cores):
  run_name = "concatenation-naive"
  run_dir = fam.get_run_dir(datadir, subst_model, "concatenation")
  shutil.rmtree(run_dir, True)
  os.makedirs(run_dir)
  supermatrix_path = os.path.join(run_dir, "supermatrix.fasta")
  partition_path = os.path.join(run_dir, "supermatrix.part")
  sites = build_supermatrix(datadir, subst_model, supermatrix_path, partition_path)
  cores = min(cores, int(sites / 500))
  start = time.time()
  
  run_raxml(subst_model, cores, run_dir, supermatrix_path, partition_path)
  time1 = (time.time() - start)
  saved_metrics.save_metrics(datadir, fam.get_run_name(run_name, subst_model), time1, "runtimes") 
  raxml_tree = os.path.join(run_dir, "concatenation.raxml.bestTree")
  dest = fam.get_species_tree(datadir, subst_model, run_name)
  shutil.copy(raxml_tree, dest)

if __name__ == "__main__":
  if (len(sys.argv) != 4):
    print("syntax: python run_concatenation.py datadir subst_model cores")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])
  run_concatenation(datadir, subst_model, cores)
  species_analyze.analyze(datadir)

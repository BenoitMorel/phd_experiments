import os
import sys
import run_raxml_supportvalues as raxml
import run_treerecs as treerecs
import run_phyldog_light as phyldog
import run_notung as notung
import run_ALE
import fam

def run_reference_methods(dataset_dir, is_dna, starting_trees, bs_trees, cores):
  print("Run pargenes light...")
  sys.stdout.flush()
  raxml.run_pargenes_and_extract_trees(dataset_dir, is_dna, 1, 0, cores, "RAxML-light", False)
  sys.stdout.flush()
  print("Run pargenes and extract trees...")
  sys.stdout.flush()
  raxml.run_pargenes_and_extract_trees(dataset_dir, is_dna, starting_trees, bs_trees, cores)
  sys.stdout.flush()
  print("Run treerecs...")
  sys.stdout.flush()
  treerecs.run_treerecs_on_families(dataset_dir, is_dna, cores)
  sys.stdout.flush()
  print("Run phyldog...")
  sys.stdout.flush()
  phyldog.run_phyldog_on_families(dataset_dir, is_dna, cores)
  print("Run notung...")
  sys.stdout.flush()
  threshold = 80
  notung.run_notung_on_families(dataset_dir, threshold, cores)
  print("Run ALE...")
  sys.stdout.flush()
  run_ALE.run_exabayes_and_ALE(dataset_dir, is_dna, cores)
  sys.stdout.flush()


if __name__ == "__main__":
  if (len(sys.argv) != 6):
    print("syntax: python run_raxml_all.py dataset_dir is_dna starting_trees bs_trees cores")
    sys.exit(1)
  dataset_dir = sys.argv[1]
  is_dna = (int(sys.argv[2]) != 0)
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
  run_reference_methods(dataset_dir, is_dna, starting_trees, bs_trees, cores)

import os
import sys
import run_raxml_supportvalues as raxml
import run_treerecs as treerecs
import run_phyldog_light as phyldog
import run_notung as notung
import run_stag

if __name__ == "__main__":
  if (len(sys.argv) != 6):
    print("syntax: python run_raxml_all.py dataset_dir is_dna starting_trees bs_trees cores")
    sys.exit(1)
  dataset_dir = sys.argv[1]
  is_dna = (int(sys.argv[2]) != 0)
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
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
  phyldog.run_phyldog_light_on_families(dataset_dir, is_dna, cores)
  print("Run notung...")
  sys.stdout.flush()
  threshold = 50
  notung.run_notung_on_families(dataset_dir, threshold, cores)


import os
import sys
import run_raxml_supportvalues as raxml
import run_treerecs as treerecs
import run_phyldog_light as phyldog
import run_notung as notung
import run_ALE
import fam
import run_generax
import analyze_dataset

class RunFilter():
  def __init__(self, raxml = True, pargenes = True, treerecs = True, phyldog = True, notung = True, generax = True, ALE = True, analyze = True):
    self.raxml = raxml
    self.pargenes = pargenes
    self.treerecs = treerecs
    self.phyldog = phyldog
    self.notung = notung
    self.generax = generax
    self.ALE = ALE
    self.analyze = analyze
  
  def disable_all(self):
    self.raxml = False
    self.pargenes = False
    self.treerecs = False
    self.phyldog = False
    self.notung = False
    self.generax = False
    self.ALE = False  
    self.analyze = False

def run_reference_methods(dataset_dir, is_dna, starting_trees, bs_trees, cores, run_filter = RunFilter()):
  if (run_filter.raxml):
    print("Run pargenes light...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(dataset_dir, is_dna, 1, 0, cores, "RAxML-light", False)
    sys.stdout.flush()
  if (run_filter.pargenes):
    print("Run pargenes and extract trees...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(dataset_dir, is_dna, starting_trees, bs_trees, cores)
    sys.stdout.flush()
  if (run_filter.treerecs):
    print("Run treerecs...")
    sys.stdout.flush()
    try:
      treerecs.run_treerecs_on_families(dataset_dir, is_dna, cores)
    except Exception as exc:
      print("Failed running Treerecs")
      print(exc)
    sys.stdout.flush()
  if (run_filter.phyldog):
    print("Run phyldog...")
    sys.stdout.flush()
    try:
      phyldog.run_phyldog_on_families(dataset_dir, is_dna, cores)
    except Exception as exc:
      print("Failed running Phyldog")
      print(exc)
  if (run_filter.notung):
    print("Run notung...")
    sys.stdout.flush()
    threshold = 80
    try:
      notung.run_notung_on_families(dataset_dir, threshold, cores)
    except Exception as exc:
      print("Failed running Notung")
      print(exc)
    threshold = 101
    try:
      notung.run_notung_on_families(dataset_dir, threshold, cores)
    except Exception as exc:
      print("Failed running Notung")
      print(exc)

  if (run_filter.generax):
    print("Run Generax")
    sys.stdout.flush()
    try:
      run_generax.run_generax_on_families(dataset_dir, is_dna, cores)
    except Exception as exc:
      print("Failed running GeneRax")
      print(exc)
    sys.stdout.flush()
  if (run_filter.ALE):
    print("Run ALE...")
    sys.stdout.flush()
    try:
      run_ALE.run_exabayes_and_ALE(dataset_dir, is_dna, cores)
    except Exception as exc:
      print("Failed running ALE")
      print(exc)
    sys.stdout.flush()
  if (run_filter.analyze):
    print("Run analyze...")
    sys.stdout.flush()
    try:
      analyze_dataset.analyze(dataset_dir)
    except Exception as exc:
      print("Failed running analyze")
      print(exc)

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

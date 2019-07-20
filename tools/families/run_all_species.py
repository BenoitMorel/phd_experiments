import os
import sys
import run_raxml_supportvalues as raxml
import fam
import run_generax
import run_stag
import species_analyze
import run_speciesrax
import run_phyldog

class SpeciesRunFilter():
  
  def __init__(self, pargenes = True, stag = True, phyldog = True, speciesrax = True, analyze = True):
    self.pargenes = pargenes
    self.stag = stag
    self.phyldog = phyldog
    self.speciesrax = speciesrax
    
    self.analyze = analyze


  def disable_all(self):
    self.pargenes = False
    self.stag = False
    self.phyldog = False
    self.speciesrax = False
    self.analyze = False

def with_transfers(dataset_dir):
  return float(dataset_dir.split("_")[-2][1:]) != 0.0

def run_reference_methods(dataset_dir, subst_model, cores, run_filter = SpeciesRunFilter()):
  if (run_filter.pargenes):
    print("Run pargenes...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(dataset_dir, subst_model, 20, 10, cores, "pargenes", True)
    sys.stdout.flush()
  if (run_filter.stag):
    print("Run Stag")
    try:
      run_stag.run_stag(dataset_dir, subst_model)
    except Exception as exc:
      print("Failed running STAG")
      print(exc)
  if (run_filter.phyldog):
    print("Run Phyldgo")
    try:
      run_phyldog.run_phyldog_on_families(dataset_dir, subst_model, cores, True)
    except Exception as exc:
      print("Failed running Phyldog")
      print(exc)
  if (run_filter.speciesrax):
    print("Run SpeciesRax")
    try:
      dtl = with_transfers(dataset_dir)
      run_speciesrax.run_speciesrax_on_families(dataset_dir, subst_model, cores, dl = (not dtl), dtl = dtl)
    except Exception as exc:
      print("Failed running speciesrax")
      print(exc)

  if (run_filter.analyze):
    print("Run analyze...")
    sys.stdout.flush()
    try:
      species_analyze.analyze(dataset_dir)
    except Exception as exc:
      print("Failed running analyze")
      print(exc)

if __name__ == "__main__":
  if (len(sys.argv) != 6):
    print("syntax: python run_raxml_all.py dataset_dir subst_model starting_trees bs_trees cores")
    sys.exit(1)
  dataset_dir = sys.argv[1]
  subst_model = sys.argv[2]
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
  run_reference_methods(dataset_dir, subst_model, starting_trees, bs_trees, cores)


import os
import sys
import run_raxml_supportvalues as raxml
import fam
import run_generax
import run_stag
import species_analyze
import run_speciesrax
import run_phyldog
import run_duptree
import run_guenomu
import run_astral_multi
import run_fastrfs

def printFlush(msg):
  print(msg)
  sys.stdout.flush()

class SpeciesRunFilter():
  
  def __init__(self):
    self.pargenes = True
    self.stag = True
    self.duptree = True
    self.fastrfs = True
    self.astral = True
    self.speciesraxfastdl = True
    self.speciesraxfastdtl = True
    self.speciesraxslowdl = True
    self.speciesraxslowdtl = True
    self.phyldog = True
    self.guenomu = False
    self.analyze = True
    

  def disable_all(self):
    self.pargenes = False
    self.stag = False
    self.fastrfs = False
    self.phyldog = False
    self.duptree = False
    self.astral = False
    self.speciesraxfastdl = False
    self.speciesraxfastdtl = False
    self.speciesraxslowdl = False
    self.speciesraxslowdtl = False
    self.guenomu = False
    #self.analyze = False
  
  def enable_fast_methods(self):
    self.disable_all()
    self.pargenes = True
    self.stag = True
    self.duptree = True
    self.fastrfs = True
    self.astral = True
    self.speciesraxfastdl = True
    self.speciesraxfastdtl = True


def with_transfers(dataset_dir):
  return float(dataset_dir.split("_")[-2][1:]) != 0.0

def run_reference_methods(dataset_dir, subst_model, cores, run_filter = SpeciesRunFilter()):
  if (run_filter.pargenes):
    printFlush("Run pargenes...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(dataset_dir, subst_model, 2, 2, cores, "pargenes", True)
    sys.stdout.flush()
  if (run_filter.stag):
    printFlush("Run Stag")
    try:
      run_stag.run_stag(dataset_dir, subst_model)
    except Exception as exc:
      printFlush("Failed running STAG\n" + str(exc))
  if (run_filter.duptree):
    printFlush("Run Duptree")
    try:
      run_duptree.run_duptree(dataset_dir, subst_model)
    except Exception as exc:
      printFlush("Failed running DupTree\n" + str(exc))
  if (run_filter.fastrfs):
    printFlush("Run FastRFS")
    try:
      run_fastrfs.run_fastrfs(dataset_dir, subst_model)
    except Exception as exc:
      printFlush("Failed running FastRFS\n" + str(exc))
  if (run_filter.astral):
    printFlush("Run Astral")
    try:
      run_astral_multi.run_astral(dataset_dir, subst_model)
    except Exception as exc:
      printFlush("Failed running Astral\n" + str(exc))
  if (run_filter.speciesraxfastdl or run_filter.speciesraxfastdtl):
    printFlush("Run SpeciesRaxFast")
    try:
      run_speciesrax.run_speciesrax_on_families(dataset_dir, subst_model, cores, dl = run_filter.speciesraxfastdl, dtl = run_filter.speciesraxfastdtl, slow = False, strategy = "SPR")
      run_speciesrax.run_speciesrax_on_families(dataset_dir, subst_model, cores, dl = False, dtl = run_filter.speciesraxfastdtl, slow = False, strategy = "TRANSFERS")
      run_speciesrax.run_speciesrax_on_families(dataset_dir, subst_model, cores, dl = False, dtl = run_filter.speciesraxfastdtl, slow = False, strategy = "HYBRID")
    except Exception as exc:
      printFlush("Failed running speciesrax\n" + str(exc))
  if (run_filter.speciesraxslowdl or run_filter.speciesraxslowdtl):
    printFlush("Run SpeciesRaxSlow")
    try:
      run_speciesrax.run_speciesrax_on_families(dataset_dir, subst_model, cores, dl = run_filter.speciesraxslowdl, dtl = run_filter.speciesraxslowdtl, slow = True)
    except Exception as exc:
      printFlush("Failed running speciesrax\n" + str(exc))
  if (run_filter.phyldog):
    printFlush("Run Phyldgo")
    try:
      run_phyldog.run_phyldog_on_families(dataset_dir, subst_model, cores, True)
    except Exception as exc:
      printFlush("Failed running Phyldog\n" + str(exc))
  if (run_filter.guenomu):
    printFlush("Run Guenomu")
    try:
      run_guenomu.run_guenomu(dataset_dir, subst_model, cores)
    except Exception as exc:
      printFlush("Failed running Guenomu\n" + str(exc))

  if (run_filter.analyze):
    printFlush("Run analyze...")
    sys.stdout.flush()
    try:
      species_analyze.analyze(dataset_dir)
    except Exception as exc:
      printFlush("Failed running analyze\n" + str(exc))

if __name__ == "__main__":
  if (len(sys.argv) != 6):
    printFlush("syntax: python run_raxml_all.py dataset_dir subst_model starting_trees bs_trees cores")
    sys.exit(1)
  dataset_dir = sys.argv[1]
  subst_model = sys.argv[2]
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
  run_reference_methods(dataset_dir, subst_model, starting_trees, bs_trees, cores)


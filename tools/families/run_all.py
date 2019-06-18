import os
import sys
import run_raxml_supportvalues as raxml
import run_treerecs as treerecs
import run_phyldog as phyldog
import run_notung as notung
import run_ALE
import fam
import run_generax
import eval_generax_likelihood
import rf_cells

class RunFilter():
  
  def __init__(self, raxml = True, pargenes = True, treerecs = True, phyldog = True, notung = True, generax = True, ALE = True, eval_joint_ll = True, analyze = True):
    self.raxml = raxml
    self.pargenes = pargenes
    self.treerecs = treerecs
    self.phyldog = phyldog
    self.notung = notung
    self.generax = generax
    self.ALE = ALE
    self.eval_joint_ll = eval_joint_ll
    self.analyze = analyze
    self.EXA_runs = -1 # use default values in run_ALE.py 
    self.EXA_chains = -1 # use default values in run_ALE.py 
    self.EXA_frequencies = -1
    self.EXA_generations = -1
    self.EXA_burnin = -1

  def disable_all(self):
    self.raxml = False
    self.pargenes = False
    self.treerecs = False
    self.phyldog = False
    self.notung = False
    self.generax = False
    self.ALE = False  
    self.eval_joint_ll = False
    self.analyze = False

def run_reference_methods(dataset_dir, subst_model, starting_trees, bs_trees, cores, run_filter = RunFilter()):
  if (run_filter.raxml):
    print("Run pargenes light...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(dataset_dir, subst_model, 1, 0, cores, "RAxML-light", False)
    sys.stdout.flush()
  if (run_filter.pargenes):
    print("Run pargenes and extract trees...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(dataset_dir, subst_model, starting_trees, bs_trees, cores)
    sys.stdout.flush()
  if (run_filter.treerecs):
    print("Run treerecs...")
    sys.stdout.flush()
    try:
      treerecs.run_treerecs_on_families(dataset_dir, subst_model, cores)
    except Exception as exc:
      print("Failed running Treerecs")
      print(exc)
    sys.stdout.flush()
  if (run_filter.phyldog):
    print("Run phyldog...")
    sys.stdout.flush()
    try:
      phyldog.run_phyldog_on_families(dataset_dir, subst_model, cores)
    except Exception as exc:
      print("Failed running Phyldog")
      print(exc)
  if (run_filter.notung):
    print("Run notung...")
    sys.stdout.flush()
    threshold = 80
    try:
      notung.run_notung_on_families(dataset_dir, subst_model,  threshold, cores)
    except Exception as exc:
      print("Failed running Notung")
      print(exc)
    #threshold = 101
    #try:
    #  notung.run_notung_on_families(dataset_dir, threshold, cores)
    #except Exception as exc:
    #  print("Failed running Notung")
    #  print(exc)

  if (run_filter.generax):
    print("Run Generax")
    sys.stdout.flush()
    try:
      run_generax.run_generax_on_families(dataset_dir, subst_model, cores)
    except Exception as exc:
      print("Failed running GeneRax")
      print(exc)
    sys.stdout.flush()
  if (run_filter.ALE):
    print("Run ALE...")
    sys.stdout.flush()
    try:
      run_ALE.run_mrbayes_and_ALE(dataset_dir, subst_model, cores, chains = run_filter.EXA_chains, runs = run_filter.EXA_runs, frequency = run_filter.EXA_frequencies, generations = run_filter.EXA_generations, burnin = run_filter.EXA_burnin)
    except Exception as exc:
      print("Failed running ALE")
      print(exc)
    sys.stdout.flush()
  if (run_filter.eval_joint_ll):
    print("Evaluating joint likelihoods...")
    try:
      eval_generax_likelihood.eval_and_save_likelihood(dataset_dir, "all", False, subst_model, cores)  
      eval_generax_likelihood.eval_and_save_likelihood(dataset_dir, "all", True, subst_model, cores)  
    except Exception as exc:
      print("Failed evaluating joint likelihoods")
      print(exc)
    sys.stdout.flush()
  if (run_filter.analyze):
    print("Run analyze...")
    sys.stdout.flush()
    try:
      rf_cells.analyze(dataset_dir)
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

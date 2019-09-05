import os
import sys
import run_raxml_supportvalues as raxml
import run_treerecs as treerecs
import run_phyldog as phyldog
import run_notung as notung
import run_eccetera
import run_deleterious
import run_ALE
import fam
import run_generax
import eval_generax_likelihood
import rf_cells
import run_mrbayes

class RunFilter():
  
  def __init__(self, raxml = True, pargenes = True, treerecs = True, phyldog = True, notung = True, eccetera = True, deleterious = True, generax = True, mrbayes = True, ALE = True, eval_joint_ll = True, analyze = True):
    self.raxml = raxml
    self.pargenes = pargenes
    self.treerecs = treerecs
    self.phyldog = phyldog
    self.notung = notung
    self.eccetera = eccetera
    self.deleterious = deleterious
    self.generax = generax
    self.generaxrec = 0
    self.mrbayes = mrbayes
    self.ALE = ALE
    self.eval_joint_ll = eval_joint_ll
    self.analyze = analyze
    self.mb_runs = 2  
    self.mb_chains = 4 
    self.mb_frequencies = 1000
    self.mb_generations = 1000000
    self.mb_burnin = 100
    self.debug = False
    self.rm_mrbayes = True

  def disable_all(self):
    self.raxml = False
    self.pargenes = False
    self.treerecs = False
    self.phyldog = False
    self.notung = False
    self.eccetera = False
    self.deleterious = False
    self.generax = False
    self.generaxrec = 0
    self.ALE = False  
    self.eval_joint_ll = False
    self.analyze = False
    self.mrbayes = False

def run_reference_methods(datadir, subst_model, starting_trees, bs_trees, cores, run_filter = RunFilter()):
  if (run_filter.raxml):
    print("Run pargenes light...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(datadir, subst_model, 1, 0, cores, "RAxML-light", False)
    sys.stdout.flush()
  if (run_filter.pargenes):
    print("Run pargenes and extract trees...")
    sys.stdout.flush()
    raxml.run_pargenes_and_extract_trees(datadir, subst_model, starting_trees, bs_trees, cores)
    sys.stdout.flush()
  if (run_filter.treerecs):
    print("Run treerecs...")
    sys.stdout.flush()
    try:
      treerecs.run_treerecs_on_families(datadir, subst_model, cores)
    except Exception as exc:
      print("Failed running Treerecs")
      print(exc)
    sys.stdout.flush()
  if (run_filter.phyldog):
    print("Run phyldog...")
    sys.stdout.flush()
    try:
      phyldog.run_phyldog_on_families(datadir, subst_model, cores)
    except Exception as exc:
      print("Failed running Phyldog")
      print(exc)
  if (run_filter.notung):
    print("Run notung...")
    sys.stdout.flush()
    threshold = 90
    try:
      notung.run_notung_on_families(datadir, subst_model,  threshold, cores)
    except Exception as exc:
      print("Failed running Notung")
      print(exc)

  if (run_filter.generax):
    print("Run Generax")
    sys.stdout.flush()
    try:
      run_generax.run_generax_on_families(datadir, subst_model, cores)
    except Exception as exc:
      print("Failed running GeneRax")
      print(exc)
    sys.stdout.flush()
  if (run_filter.generaxrec > 0):
    print("Run Generax Rec")
    sys.stdout.flush()
    try:
      run_generax.run_generax_on_families(datadir, subst_model, cores, True, True, True, True, run_filter.generaxrec)
    except Exception as exc:
      print("Failed running GeneRax")
      print(exc)
    sys.stdout.flush()
  if (run_filter.mrbayes):
    print("Run mrbayes...")
    sys.stdout.flush()
    try:
      run_mrbayes.run_mrbayes_on_families(datadir, run_filter.mb_generations, run_filter.mb_frequencies, run_filter.mb_runs, run_filter.mb_chains, run_filter.mb_burnin, subst_model, cores)
    except Exception as exc:
      print("Failed running mrbayes")
      print(exc)
  if (run_filter.ALE):
    print("Run ALE...")
    sys.stdout.flush()
    try:
      run_ALE.run_ALE(datadir, subst_model, cores)
    except Exception as exc:
      print("Failed running ALE")
      print(exc)
    sys.stdout.flush()
  if (run_filter.eccetera):
    print("Run eccetera...")
    sys.stdout.flush()
    threshold = 70
    try:
      run_eccetera.run_eccetera_on_families(datadir, subst_model,  threshold, cores)
    except Exception as exc:
      print("Failed running Eccetera")
      print(exc)
  if (run_filter.deleterious):
    print("Run deleterious...")
    sys.stdout.flush()
    try:
      run_deleterious.run_deleterious_on_families(datadir, subst_model, cores)
    except Exception as exc:
      print("Failed running Deleterious")
      print(exc)
  if (run_filter.eval_joint_ll):
    print("Evaluating joint likelihoods...")
    try:
      eval_generax_likelihood.eval_and_save_likelihood(datadir, "all", False, subst_model, cores)  
      eval_generax_likelihood.eval_and_save_likelihood(datadir, "all", True, subst_model, cores)  
    except Exception as exc:
      print("Failed evaluating joint likelihoods")
      print(exc)
    sys.stdout.flush()
  if (run_filter.analyze):
    print("Run analyze...")
    sys.stdout.flush()
    try:
      rf_cells.analyze(datadir)
    except Exception as exc:
      print("Failed running analyze")
      print(exc)
  if (run_filter.rm_mrbayes):
    try:
      print("Removing mrbayes files...")
      run_mrbayes.remove_mrbayes_run(datadir, subst_model)
    except:
      pass

if __name__ == "__main__":
  if (len(sys.argv) != 6):
    print("syntax: python run_raxml_all.py datadir subst_model starting_trees bs_trees cores")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
  run_reference_methods(datadir, subst_model, starting_trees, bs_trees, cores)

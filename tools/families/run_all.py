import os
import sys
import run_raxml_supportvalues as raxml
import run_treerecs as treerecs
import run_treefixdtl as treefix
import run_phyldog as phyldog
import run_notung as notung
import run_eccetera
import run_deleterious
import run_ALE
import fam
import run_generax
import eval_generax_likelihood
import fast_rf_cells
import run_mrbayes
import run_stag
import species_analyze
import run_speciesrax

def printFlush(msg):
  print(msg)
  sys.stdout.flush()


class RunFilter():
  def __init__(self, gene_inference = True, species_inference = False):
    self.disable_all()
    if (gene_inference):
      self.raxml = True
      self.pargenes = True
      self.treerecs = True
      self.treefix = False
      self.phyldog = True
      self.notung = True
      self.eccetera = True
      self.generax = True
      self.mrbayes = True
      self.ALE = True
      self.eval_joint_ll = True
      self.analyze = True
      self.rm_mrbayes = True
    if (species_inference):
      self.raxml = True 
      self.pargenes = True
      self.stag = True
      self.phyldog_species = True
      self.speciesrax = True
      self.analyze_species = True

    self.pargenes_starting_trees = 20
    self.pargenes_bootstrap_trees = 100
    self.mb_runs = 2  
    self.mb_chains = 4 
    self.mb_frequencies = 1000
    self.mb_generations = 1000000
    self.mb_burnin = 100
    self.deleterious = False
    self.debug = False
    self.dated_ALE = False

  def disable_all(self):
    self.raxml = False
    self.pargenes = False
    self.treerecs = False
    self.treefix = False
    self.phyldog = False
    self.notung = False
    self.eccetera = False
    self.deleterious = False
    self.generax = False
    self.ALE = False  
    self.eval_joint_ll = False
    self.analyze = False
    self.mrbayes = False
    self.dated_ALE = False
    self.stag = False
    self.phyldog_species = False
    self.speciesrax = False
    self.analyze_species = False

  def run_reference_methods(self, datadir, subst_model, cores):
    print("*************************************")
    print("Run tested gene tree inference tools for dataset " + datadir)
    print("*************************************")
    if (len(datadir.split("/")) == 1):
      datadir = fam.get_datadir(datadir) 
    save_stdout = sys.stdout
    redirected_file = os.path.join(datadir, "runs", "logs_run_all_genes." + subst_model + ".txt")
    print("Redirected logs to " + redirected_file)
    sys.stdout.flush()
    #sys.stdout = open(redirected_file, "w")
    if (self.raxml):
      printFlush("Run pargenes light...")
      raxml.run_pargenes_and_extract_trees(datadir, subst_model, 1, 0, cores, "RAxML-light", False)
    
    if (self.pargenes):
      printFlush("Run pargenes and extract trees...")
      raxml.run_pargenes_and_extract_trees(datadir, subst_model, self.pargenes_starting_trees, self.pargenes_bootstrap_trees, cores)
    
    if (self.treerecs):
      printFlush("Run treerecs...")
      try:
        treerecs.run_treerecs_on_families(datadir, subst_model, cores)
      except Exception as exc:
        printFlush("Failed running Treerecs\n" + str(exc))
    if (self.treefix):
      printFlush("Run treefix...")
      try:
        treefix.run_treefix_on_families(datadir, subst_model, cores)
      except Exception as exc:
        printFlush("Failed running TreeFixs\n" + str(exc))
    
    if (self.phyldog):
      printFlush("Run phyldog...")
      try:
        phyldog.run_phyldog_on_families(datadir, subst_model, cores)
      except Exception as exc:
        printFlush("Failed running Phyldog\n" + str(exc))
    
    if (self.notung):
      printFlush("Run notung...")
      threshold = 90
      try:
        notung.run_notung_on_families(datadir, subst_model,  threshold, cores)
      except Exception as exc:
        printFlush("Failed running Notung\n" + str(exc))

    if (self.generax):
      printFlush("Run Generax")
      try:
        run_generax.run_generax_on_families(datadir, subst_model, cores)
      except Exception as exc:
        printFlush("Failed running GeneRax\n" + str(exc))
    
    if (self.stag):
      printFlush("Run Stag")
      try:
        run_stag.run_stag(datadir, subst_model)
      except Exception as exc:
        printFlush("Failed running STAG\n" + str(exc))
    
    
    if (self.speciesrax):
      printFlush("Run SpeciesRax")
      try:
        run_generax.run_generax_on_families(datadir, subst_model, cores, raxml = False, random = True, optimize_species = True)
      except Exception as exc:
        printFlush("Failed running speciesrax\n" + str(exc))
    
    if (self.phyldog_species):
      printFlush("Run Phyldog species")
      try:
        phyldog.run_phyldog_on_families(datadir, subst_model, cores, True)
      except Exception as exc:
        printFlush("Failed running Phyldog species\n" + str(exc))


    if (self.mrbayes):
      printFlush("Run mrbayes...")
      try:
        instance = run_mrbayes.MrbayesInstance(datadir, subst_model, self.mb_runs, self.mb_chains, self.mb_generations, self.mb_frequencies, self.mb_burnin) 
        run_mrbayes.run_mrbayes_on_families(instance, cores, False)
      except Exception as exc:
        printFlush("Failed running mrbayes\n" + str(exc))
    
    if (self.ALE):
      printFlush("Run ALE...")
      try:
        run_ALE.run_ALE(datadir, subst_model, cores)
      except Exception as exc:
        printFlush("Failed running ALE\n" + str(exc))
    if (self.dated_ALE):
      printFlush("Run dated ALE...")
      try:
        run_ALE.run_ALE(datadir, subst_model, cores, True)
      except Exception as exc:
        printFlush("Failed running dated ALE\n" + str(exc))
    if (self.eccetera):
      printFlush("Run eccetera...")
      threshold = 70
      try:
        run_eccetera.run_eccetera_on_families(datadir, subst_model,  threshold, cores)
      except Exception as exc:
        printFlush("Failed running Eccetera\n" + str(exc))
    
    if (self.deleterious):
      printFlush("Run deleterious...")
      try:
        run_deleterious.run_deleterious_on_families(datadir, subst_model, cores)
      except Exception as exc:
        printFlush("Failed running Deleterious\n" + str(exc))
    
    if (self.eval_joint_ll):
      printFlush("Evaluating joint likelihoods...")
      try:
        eval_generax_likelihood.eval_and_save_likelihood(datadir, "all", False, subst_model, cores)  
        eval_generax_likelihood.eval_and_save_likelihood(datadir, "all", True, subst_model, cores)  
      except Exception as exc:
        printFlush("Failed evaluating joint likelihoods\n" + str(exc))
    
    if (self.analyze):
      printFlush("Run analyze...")
      try:
        fast_rf_cells.analyze(datadir, "all", cores)
      except Exception as exc:
        printFlush("Failed running analyze\n" + str(exc))
    
    if (self.analyze_species):
      printFlush("Run analyze species...")
      try:
        species_analyze.analyze(datadir)
      except Exception as exc:
        printFlush("Failed running analyze\n" + str(exc))
    
    if (self.rm_mrbayes):
      try:
        printFlush("Removing mrbayes files...")
        run_mrbayes.remove_mrbayes_run(datadir, subst_model)
      except:
        printFlush("Failed removing mrbayes files\n", str(exc))
    print("Done")
    sys.stdout = save_stdout
    
  def run_all_reference_methods(self, datasets, subst_model, cores = 40):
    for dataset in datasets:
       self.run_reference_methods(dataset, subst_model, cores)

if __name__ == "__main__":
  if (len(sys.argv) != 6):
    printFlush("syntax: python run_raxml_all.py datadir subst_model starting_trees bs_trees cores")
    sys.exit(1)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  starting_trees = sys.argv[3]
  bs_trees = sys.argv[4]
  cores = int(sys.argv[5])
  run_reference_methods(datadir, subst_model, starting_trees, bs_trees, cores)

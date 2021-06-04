import os
import sys
import pickle
import uuid
import tempfile
import time
import run_raxml_supportvalues as raxml
import fam
import run_generax
import run_stag
import species_analyze
import run_speciesrax
import run_tegrator
import run_phyldog
import run_duptree
import run_guenomu
import run_astrid
import run_astral_multi
import run_astral_pro
import run_fastrfs
import run_fastmulrfs
import run_njst  
import run_njrax
import run_concatenation
import run_orthogenerax
import run_generax_selector
import run_mrbayes
import run_fasttree
import shutil
sys.path.insert(0, os.path.join("scripts"))
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import fam_data


def printFlush(msg):
  print(msg)
  sys.stdout.flush()

class SpeciesRunFilter():
  
  def __init__(self):
    self.generate = False
    self.pargenes = True
    self.pargenes_starting_trees = 1
    self.pargenes_bootstrap_trees = 0
    self.fasttree = False
    self.mrbayes = False
    self.mb_runs = 1  
    self.mb_chains = 4 
    self.mb_frequencies = 1000
    self.mb_generations = 10000
    self.mb_burnin = 1
    self.starting_gene_trees = ["raxml-ng"]
    self.orthogenerax = True
    self.concatenation_min = True
    self.concatenation_max = True
    self.stag = True
    self.duptree = True
    self.fastrfs = True
    self.fastmulrfs = True
    self.njrax = True    
    self.cherry = True    
    self.cherrypro = True    
    self.astrid = True    
    self.astral = True
    self.generaxselect = True
    self.generaxselectfam = True
    self.astralpro = True
    self.speciesrax = True
    self.speciesraxparsidl = True
    self.speciesraxprune = True
    self.speciesraxnofam = True
    self.speciesraxperfamily = True
    self.speciesraxbench = True
    self.genetegratorbench = True
    self.phyldog = True
    self.guenomu = False
    self.analyze = True
    self.cleanup = False
    self.verbose = False

  def disable_all(self):
    self.pargenes = False
    self.fasttree = False
    self.mrbayes = False
    self.orthogenerax = False
    self.concatenation_min = False
    self.concatenation_max = False
    self.stag = False
    self.fastrfs = False
    self.fastmulrfs = False
    self.phyldog = False
    self.duptree = False
    self.njrax = False    
    self.cherry = False    
    self.cherrypro = False    
    self.njst = False    
    self.astrid = False
    self.astral = False
    self.generaxselect = False
    self.generaxselectfam = False
    self.astralpro = False
    self.speciesrax = False
    self.speciesraxparsidl = False
    self.speciesraxnofam = False
    self.speciesraxprune = False
    self.speciesraxperfamily = False
    self.genetegratorbench = False
    self.speciesraxbench = False
    self.guenomu = False
    self.cleanup = False
    #self.analyze = False
  
  def enable_fast_methods(self):
    self.disable_all()
    self.pargenes = True
    self.fasttree = True
    self.orthogenerax = True
    self.concatenation_min = True
    self.concatenation_max = True
    self.stag = True
    self.duptree = True
    self.fastrfs = True
    self.fastmulrfs = True
    self.njrax = True   
    self.cherry = True   
    self.cherrypro = True   
    self.njst = True   
    self.astrid = True
    self.speciesrax = True

  def launch_reference_methods(self, datadir, subst_model, cores, launch_mode):
    misc_dir = fam.get_misc_dir(datadir)
    run_filter_file_path = os.path.join(misc_dir, str(uuid.uuid4()))
    pickle.dump(self, open(run_filter_file_path, "wb"))
    command = []
    command.append(exp.python())
    command.append(os.path.realpath(__file__))
    command.append(datadir)
    command.append(subst_model)
    command.append(str(cores))
    command.append(run_filter_file_path)
    submit_path = os.path.join(misc_dir, "submit_" + str(uuid.uuid4()) + ".sh")
    print("Submit path : " + submit_path)
    print("Saving pickle in " + run_filter_file_path)
    exp.submit(submit_path, " ".join(command), cores, launch_mode) 
    


  def run_reference_methods(self, datadir, subst_model, cores, launch_mode = "normal"):
    if (self.generate):
      if (not os.path.isdir(datadir)):
        print("Generating " + datadir)
        dataset = os.path.basename(datadir)
        fam_data.generate_dataset(dataset)

    if (launch_mode != "normal"):
      self.launch_reference_methods(datadir, subst_model, cores, launch_mode)
      return
    print("*************************************")
    print("Run tested species inference tools for dataset " + datadir)
    print("*************************************")
    if (len(datadir.split("/")) == 1):
      datadir = fam.get_datadir(datadir) 
    save_sdtout = sys.stdout
    try:
      os.makedirs(os.path.join(datadir, "runs"))
    except:
      pass
    if (not self.verbose):
      redirected_file = os.path.join(datadir, "runs", "logs_run_all_species." + subst_model + ".txt")
      print("Redirected logs to " + redirected_file)
      sys.stdout.flush()
      sys.stdout = open(redirected_file, "w")
    if (self.pargenes and subst_model != "true"):
      printFlush("Run pargenes...")
      sys.stdout.flush()
      raxml.run_pargenes_and_extract_trees(datadir, subst_model, self.pargenes_starting_trees, self.pargenes_bootstrap_trees, cores, "pargenes", True)
      sys.stdout.flush()
    if (self.fasttree and subst_model != "true"):
      printFlush("Run FastTree...")
      sys.stdout.flush()
      run_fasttree.run_fasttree_on_families(datadir, subst_model, cores)
      sys.stdout.flush()
    if (self.mrbayes):
      printFlush("Run mrbayes...")
      try:
        instance = run_mrbayes.MrbayesInstance(datadir, subst_model, self.mb_runs, self.mb_chains, self.mb_generations, self.mb_frequenciesfrequency, self.mb_burnin) 
        run_mrbayes.run_mrbayes_on_families(instance, cores, False)
      except Exception as exc:
        printFlush("Failed running mrbayes\n" + str(exc))
    if (self.stag):
      for gene_tree in self.starting_gene_trees:
        printFlush("Run Stag")
        try:
          run_stag.run_stag(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running STAG\n" + str(exc))
    if (self.duptree):
      printFlush("Run Duptree")
      for gene_tree in self.starting_gene_trees:
        try:
          run_duptree.run_duptree(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running DupTree\n" + str(exc))
    if (self.fastmulrfs):
      printFlush("Run FastMulRFS")
      for gene_tree in self.starting_gene_trees:
        try:
          run_fastmulrfs.run_fastmulrfs(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running FastMulRFS\n" + str(exc))
    if (self.fastrfs):
      printFlush("Run FastRFS")
      for gene_tree in self.starting_gene_trees:
        try:
          run_fastrfs.run_fastrfs(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running FastRFS\n" + str(exc))
    if (self.njrax):
      printFlush("Run NJrax")
      for gene_tree in self.starting_gene_trees:
        try:
          run_njrax.run_njrax(datadir, "MiniNJ", gene_tree, subst_model)
          run_njrax.run_njrax(datadir, "WMiniNJ", gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running NJrax with " + gene_tree + "\n" + str(exc))
    if (self.cherry):
      printFlush("Run Cherry")
      for gene_tree in self.starting_gene_trees:
        try:
          run_njrax.run_njrax(datadir, "Cherry", gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running Cherry with " + gene_tree + "\n" + str(exc))
    if (self.cherrypro):
      printFlush("Run CherryPro")
      for gene_tree in self.starting_gene_trees:
        try:
          run_njrax.run_njrax(datadir, "CherryPro", gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running CherryPro with " + gene_tree + "\n" + str(exc))
    if (self.njst):
      printFlush("Run NJst")
      for gene_tree in self.starting_gene_trees:
        try:
          run_njrax.run_njrax(datadir, "NJst", gene_tree, subst_model)
          run_njrax.run_njrax(datadir, "Ustar", gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running NJst with " + gene_tree + "\n" + str(exc))
    if (self.astrid):
      printFlush("Run Astrid")
      try:
        run_astrid.run_astrid(datadir, "raxml-ng", subst_model)
      except Exception as exc:
        printFlush("Failed running Astrid\n" + str(exc))
    if (self.astral):
      printFlush("Run Astral")
      try:
        run_astral_multi.run_astral(datadir, "raxml-ng", subst_model)
      except Exception as exc:
        printFlush("Failed running Astral\n" + str(exc))
    if (self.astralpro):
      printFlush("Run Astral-pro")
      for gene_tree in self.starting_gene_trees:
        try:
          run_astral_pro.run_astralpro(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running Astral-pro with " + gene_tree + "\n" + str(exc))
    if (self.concatenation_min and subst_model != "true"):
      printFlush("Run concatenation-min")
      try:
        run_concatenation.run_concatenation(datadir, subst_model,False,cores)
      except Exception as exc:
        printFlush("Failed running concatenation-min\n" + str(exc))
    if (self.concatenation_max):
      printFlush("Run concatenation-max"  and subst_model != "true")
      try:
        run_concatenation.run_concatenation(datadir, subst_model,True,cores)
      except Exception as exc:
        printFlush("Failed running concatenation-min\n" + str(exc))
    for gene_tree in self.starting_gene_trees:
      if (self.speciesrax):
        printFlush("Run SpeciesRaxFast")
        try:
          dataset = os.path.basename(datadir)
          run_speciesrax.run_speciesrax_on_families(dataset, gene_tree, subst_model, cores, transfers = True, strategy = "HYBRID", rates_per_family = False)
        except Exception as exc:
          printFlush("Failed running speciesrax\n" + str(exc))
      if (self.speciesraxparsidl):
        printFlush("Run SpeciesRaxParsi")
        try:
          dataset = os.path.basename(datadir)
          run_speciesrax.run_speciesrax_instance(dataset, gene_tree, 2, "speciesrax-parsidl", subst_model, "HYBRID", cores)
        except Exception as exc:
          printFlush("Failed running speciesrax prune\n" + str(exc))
      if (self.speciesraxnofam):
        printFlush("Run SpeciesRaxNoFam")
        try:
          dataset = os.path.basename(datadir)
          run_speciesrax.run_speciesrax_instance(dataset, gene_tree, True, "speciesrax-prune", subst_model, "HYBRID", cores, ["--prune-species-tree"])
        except Exception as exc:
          printFlush("Failed running speciesrax prune\n" + str(exc))
      if (self.speciesraxprune):
        printFlush("Run SpeciesRaxPrune")
        try:
          dataset = os.path.basename(datadir)
          run_speciesrax.run_speciesrax_instance(dataset, gene_tree, True, "speciesrax-prune", subst_model, "HYBRID", cores, ["--prune-species-tree", "--per-family-rates"])
        except Exception as exc:
          printFlush("Failed running speciesrax prune\n" + str(exc))
      if (self.genetegratorbench):
        printFlush("Run genetegrator bench")
        dataset = os.path.basename(datadir)
        run_tegrator.run_genetegrator_bench(dataset, "generax-MiniNJ-fam_raxml-ng", "mrbayes", subst_model, cores)
      if (self.speciesraxbench):
        printFlush("Run speciesRaxBench")
        dataset = os.path.basename(datadir)
        try:
          no_opt_args = [] #"--dtl-rates-opt", "NONE"]
          run_speciesrax.run_speciesrax_bench(dataset, "UndatedDTL", False, "MiniNJ", gene_tree, subst_model, cores, no_opt_args)
          run_speciesrax.run_speciesrax_bench(dataset, "UndatedDTL", False, "random", gene_tree, subst_model, cores, ["--seed", "1"] + no_opt_args)
        except Exception as exc:
          printFlush("Failed running speciesrax bench\n" + str(exc))

      if(self.speciesraxperfamily):
        printFlush("Run SpeciesRaxPerFamily")
        try:
          dataset = os.path.basename(datadir)
          run_speciesrax.run_speciesrax_on_families(datadir, gene_tree, subst_model, cores, transfers = True, strategy = "HYBRID", rates_per_family = True)
        except Exception as exc:
          printFlush("Failed running speciesrax per family\n" + str(exc))

    if (self.phyldog):
      printFlush("Run Phyldgo")
      try:
        run_phyldog.run_phyldog_on_families(datadir, subst_model, cores, True)
      except Exception as exc:
        printFlush("Failed running Phyldog\n" + str(exc))
    if (self.guenomu):
      printFlush("Run Guenomu")
      try:
        run_guenomu.run_guenomu(datadir, subst_model, cores)
      except Exception as exc:
        printFlush("Failed running Guenomu\n" + str(exc))
    if (self.orthogenerax  and subst_model != "true"):
      printFlush("Run OrthoGeneRax")
      try:
        #run_orthogenerax.run_orthogenerax(datadir, subst_model, "true", False, cores)
        #run_orthogenerax.run_orthogenerax(datadir, subst_model, "njrax-NJst", False, cores)
        run_orthogenerax.run_orthogenerax(datadir, subst_model, "astralpro-raxml-ng", False, cores)
        #run_orthogenerax.run_orthogenerax(datadir, subst_model, "speciesrax-dtl-raxml-HYBRID", False, cores)
      except Exception as exc:
        printFlush("Failed running orthogenerax\n" + str(exc))

    if (self.generaxselect  and subst_model != "true"):
      printFlush("Run GeneRaxSelect")
      try:
        candidates = exp.generax_selector_candidates
        run_generax_selector.select(datadir, subst_model, candidates, True, cores)
      except Exception as exc:
        printFlush("Failed running GeneRaxSelect\n" + str(exc))
    if (self.generaxselectfam  and subst_model != "true"):
      printFlush("Run GeneRaxSelectFam")
      try:
        candidates = exp.generax_selector_candidates
        run_generax_selector.select(datadir, subst_model, candidates, True, cores)
      except Exception as exc:
        printFlush("Failed running GeneRaxSelect\n" + str(exc))
    if (self.analyze):
      printFlush("Run analyze...")
      sys.stdout.flush()
      try:
        species_analyze.analyze(datadir)
      except Exception as exc:
        printFlush("Failed running analyze\n" + str(exc))
    if (self.cleanup):
      try:
        shutil.rmtree(fam.get_run_dir(datadir, subst_model))
      except Exception as exc:
        print("Failed to remove run directory: \n" + str(exc))
    sys.stdout = save_sdtout
    print("End of run_all")
    sys.stdout.flush()
  
  def run_all_reference_methods(self, datasets, subst_model, cores = 40):
    for dataset in datasets:
       self.run_reference_methods(dataset, subst_model, cores)


if __name__ == "__main__":
  if (len(sys.argv) != 5):
    print("syntax: python run_all_species.py datadir subst_model cores run_filter_pickle")
    sys.exit(1)
  

  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  cores = int(sys.argv[3])
  run_filter_pickle = sys.argv[4]
  run_filter = pickle.load(open(run_filter_pickle, "rb"))
  run_filter.run_reference_methods(datadir, subst_model, cores)


import os
import sys
import pickle
import uuid
import tempfile
import time
import run_raxml_supportvalues as raxml
import fam
import run_generax
import run_alegenerax
import run_ALE
import run_stag
import species_analyze
import grf_species_analyze
import run_speciesrax
import run_tegrator
import run_phyldog
import run_duptree
import run_stride
import run_guenomu
import run_asteroid
import run_astrid
import run_astrid_single
import run_astral_multi
import run_astral
import run_aster
import run_mp_astral
import run_astral_pro
import run_aster
import run_stells
import run_fastrfs
import run_fastmulrfs
import run_njst  
import run_njrax
import run_concatenation
import run_mrbayes
import run_fasttree
import run_dicotree
import fast_rf_cells
import shutil
import subprocess
import run_bootstrap_trees
import contract_branches
import numpy as np

sys.path.insert(0, os.path.join("scripts"))
sys.path.insert(0, os.path.join("tools", "families"))
import experiments as exp
import fam_data


def printFlush(msg):
  print(msg)
  sys.stdout.flush()

def floatstr(my_float):
  return np.format_float_positional(my_float, trim='-')

class SpeciesRunFilter():
  
  def __init__(self):
    self.generate = False
    self.disable_all()
    self.starting_gene_trees = ["raxml-ng"]
    
    # gene trees
    self.pargenes_starting_trees = 1
    self.pargenes_bootstrap_trees = 0
    self.bootstrap_trees = 0
    self.mb_runs = 1  
    self.mb_chains = 4 
    self.mb_frequencies = 1000
    self.mb_generations = 10000
    self.mb_burnin = 1
    self.minbl = -1.0
    self.analyze_gene_trees = False
    self.cleanup = False
    self.verbose = False

  def disable_all(self):
    self.generate = False
    self.pargenes = False
    self.bootstrap_trees = 0
    self.fasttree = False
    self.dicotree = False
    self.mrbayes = False
    self.generax_undated = False
    self.generax_softdated = False
    self.ale_undated = []
    self.ale_dated = []
    self.alegenerax_undated = []
    self.alegenerax_softdated = []
    self.alegenerax_softdated_gamma = []
    self.concatenation_min = False
    self.concatenation_max = False
    self.stag = False
    self.stells = False
    self.fastrfs = False
    self.fastmulrfs = False
    self.phyldog = False
    self.duptree = False
    self.stride = False
    self.njrax = False    
    self.cherry = False    
    self.cherrypro = False    
    self.njst = False    
    self.astrid = []
    self.asteroid = []
    self.astrid_single = []
    self.aster = []
    self.astral = []
    self.astral_mp = []
    self.astralpro = []
    self.astralpro2 = []
    self.speciesrax = []
    self.genetegratorbench = []
    self.speciesraxbench = []
    self.guenomu = False
    self.cleanup = False
    self.analyze = False
    self.grf_analyze = False
    self.analyze_gene_trees = False
  
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
    printFlush("HEY " + os.path.join(datadir, "runs", subst_model))
    try:
      os.makedirs(os.path.join(datadir, "runs", subst_model))
    except:
      printFlush("Cannot build " + os.path.join(datadir, "runs", subst_model))
      pass
    printFlush("ok, built!")
    if (not self.verbose):
      redirected_file = os.path.join(datadir, "runs", "logs_run_all_species." + subst_model + ".txt")
      print("Redirected logs to " + redirected_file)
      sys.stdout.flush()
      sys.stdout = open(redirected_file, "w")
    if (self.pargenes and subst_model != "true"):
      printFlush("Run pargenes...")
      sys.stdout.flush()
      raxml.run_pargenes_and_extract_trees(datadir, subst_model, 0, self.pargenes_starting_trees, self.pargenes_bootstrap_trees, cores, "pargenes", True)
      sys.stdout.flush()
    if (self.fasttree and subst_model != "true"):
      printFlush("Run FastTree...")
      sys.stdout.flush()
      assert(os.path.isdir((os.path.join(datadir, "runs", subst_model))))
      run_fasttree.run_fasttree_on_families(datadir, subst_model, cores)
      sys.stdout.flush()
    if (self.dicotree):
      printFlush("Run DiCoTree...")
      sys.stdout.flush()
      assert(os.path.isdir((os.path.join(datadir, "runs", subst_model))))
      run_dicotree.run_dicotree_on_families(datadir, subst_model, cores)
      sys.stdout.flush()
    if (self.bootstrap_trees != 0):
      run_bootstrap_trees.run_pargenes_and_extract_trees(datadir, subst_model, self.bootstrap_trees, cores)
    if (self.mrbayes):
      printFlush("Run mrbayes...")
      try:
        instance = run_mrbayes.MrbayesInstance(datadir, subst_model, self.mb_runs, self.mb_chains, self.mb_generations, self.mb_frequencies, self.mb_burnin) 
        run_mrbayes.run_mrbayes_on_families(instance, cores, False)
        instance.remove_mrbayes_run()
      except Exception as exc:
        printFlush("Failed running mrbayes\n" + str(exc))
    if (self.generax_undated):
      printFlush("Run generax undated")
      run_generax.run(datadir, subst_model, "PARENTS", cores, [])
    if (self.generax_softdated):
      printFlush("Run generax softdated")
      run_generax.run(datadir, subst_model, "SOFTDATED", cores, [])
    for gene_tree in self.alegenerax_undated:
      printFlush("Run alegenerax Undated from  " + gene_tree )
      run_alegenerax.run(datadir, gene_tree, subst_model, "PARENTS", cores, [])
    for gene_tree in self.alegenerax_softdated:
      printFlush("Run alegenerax SoftDated from  " + gene_tree )
      run_alegenerax.run(datadir, gene_tree, subst_model, "SOFTDATED", cores, [])
    for gene_tree in self.alegenerax_softdated_gamma:
      printFlush("Run alegenerax SoftDated from  " + gene_tree )
      run_alegenerax.run(datadir, gene_tree, subst_model, "SOFTDATED", cores, ["--gamma-categories", "4"])
    for gene_tree in self.ale_undated:
      printFlush("Run ALE Undated from  " + gene_tree )
      run_ALE.run_ALE(datadir, gene_tree, subst_model, cores, dated = False)
    for gene_tree in self.ale_dated:
      printFlush("Run ALE Dated from  " + gene_tree )
      run_ALE.run_ALE(datadir, gene_tree, subst_model, cores, dated = True)

      run_generax.run(datadir, subst_model, "SOFTDATED", cores, [])
    if (self.minbl > 0.0):
      for gene_tree in self.starting_gene_trees:
        print("Coucou " + gene_tree)
        gene_tree = gene_tree.split("-minbl")[0]
        print("Coucou " + gene_tree)
        printFlush("Contracting bl < " + floatstr(self.minbl))
        contract_branches.contract(datadir, gene_tree, subst_model, self.minbl, -1.0)
    if (self.stag):
      for gene_tree in self.starting_gene_trees:
        printFlush("Run Stag")
        try:
          run_stag.run_stag(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running STAG\n" + str(exc))
    if (self.stride):
      printFlush("Running stride...")
      for gene_tree in self.starting_gene_trees:
        try:
          run_stride.run_stride(datadir, "generax-MiniNJ-fam", gene_tree, subst_model)
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
    if (self.stells):
      printFlush("Run FastRFS")
      for gene_tree in self.starting_gene_trees:
        try:
          run_stells.run_stells(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running STELLS\n" + str(exc))
    if (self.njrax):
      printFlush("Run NJrax")
      for gene_tree in self.starting_gene_trees:
        try:
          run_njrax.run_njrax(datadir, "MiniNJ", gene_tree, subst_model)
          #run_njrax.run_njrax(datadir, "WMiniNJ", gene_tree, subst_model)
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
    for gene_tree in self.asteroid:
      printFlush("Run Asteroid " + gene_tree)
      run_asteroid.run_asteroid(datadir, gene_tree, subst_model, cores, ["-r", "1"])
    for gene_tree in self.astrid:
      printFlush("Run Astrid " + gene_tree)
      for gene_tree in self.starting_gene_trees:
        try:
          run_astrid.run_astrid(datadir, gene_tree, subst_model, "fastme")
        except Exception as exc:
          printFlush("Failed running Astrid\n" + str(exc))
    for gene_tree in self.astrid_single:
      printFlush("Run Astrid")
      for gene_tree in self.starting_gene_trees:
        run_astrid_single.run_astrid(datadir, gene_tree, subst_model, "fastme")
    for gene_tree in self.astral_mp:
      printFlush("Run Astral")
      run_mp_astral.run_astral(datadir, gene_tree, subst_model)
    for gene_tree in self.astral:
      printFlush("Run Astral")
      run_astral.run_astral(datadir, gene_tree, subst_model)
    for gene_tree in self.aster:
      printFlush("Run Aster")
      run_aster.run_aster(datadir, gene_tree, subst_model, cores)
    if (self.concatenation_min and subst_model != "true"):
      printFlush("Run concatenation-min")
      try:
        run_concatenation.run_concatenation(datadir, "min", subst_model, cores)
      except Exception as exc:
        printFlush("Failed running concatenation-min\n" + str(exc))
    if (self.concatenation_max):
      printFlush("Run concatenation-max"  and subst_model != "true")
      try:
        run_concatenation.run_concatenation(datadir, "max", subst_model,cores)
      except Exception as exc:
        printFlush("Failed running concatenation-min\n" + str(exc))
    for gene_tree in self.speciesrax:
      dataset = os.path.basename(datadir)
      run_speciesrax.run_speciesrax_on_families(dataset, gene_tree, subst_model, cores, transfers = True, strategy = "HYBRID", rates_per_family = True)
    for gene_tree in self.astralpro:
      printFlush("Run Astral-pro")
      run_astral_pro.run_astralpro(datadir, gene_tree, subst_model)
    for gene_tree in self.astralpro2:
      printFlush("Run Astral-pro 2")
      run_aster.run_aster(datadir, gene_tree, subst_model, True, cores)
    for gene_tree in self.genetegratorbench:
      printFlush("Run genetegrator bench")
      run_tegrator.run_genetegrator_bench(datadir, "MiniNJ", gene_tree, subst_model, cores)
      #run_tegrator.run_genetegrator_bench(datadir, "MiniNJ", gene_tree, subst_model, cores, ["--rec-model", "UndatedDL"])
    for gene_tree in self.speciesraxbench:
      printFlush("Run speciesRaxBench")
      dataset = os.path.basename(datadir)
      run_speciesrax.run_speciesrax_bench(dataset, "UndatedDTL", False, "MiniNJ", gene_tree, subst_model, cores, [])
    if (self.duptree):
      printFlush("Run Duptree")
      for gene_tree in self.starting_gene_trees:
        try:
          run_duptree.run_duptree(datadir, gene_tree, subst_model)
        except Exception as exc:
          printFlush("Failed running DupTree\n" + str(exc))
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

    if (self.analyze):
      printFlush("Run analyze...")
      sys.stdout.flush()
      try:
        species_analyze.analyze(datadir)
      except Exception as exc:
        printFlush("Run analyze gene trees...")
    if (self.grf_analyze):
      printFlush("Run grf analyze...")
      sys.stdout.flush()
      try:
        grf_species_analyze.analyze(datadir)
      except Exception as exc:
        printFlush("Run grf analyze failed...")
    if (self.analyze_gene_trees):
      try:
        run = "all"
        fast_rf_cells.analyze(datadir, run, cores)
      except Exception as exc:
        printFlush("Failed running analyze\n" + str(exc))
    if (self.cleanup):
      try:
        print("REMOVE TREE " + fam.get_run_dir(datadir, subst_model))
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


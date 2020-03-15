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
import run_phyldog
import run_duptree
import run_guenomu
import run_astrid
import run_astral_multi
import run_astral_pro
import run_fastrfs
import run_njst  
import run_njrax
import run_concatenation
import run_orthogenerax
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
    self.orthogenerax = True
    self.concatenation_min = True
    self.concatenation_max = True
    self.stag = True
    self.duptree = True
    self.fastrfs = True
    self.njrax = True    
    self.njst = True    
    self.astrid = True    
    self.astral = True
    self.astralpro = True
    self.speciesraxfastdl = True
    self.speciesraxfastdtl = True
    self.speciesraxprune = True
    self.speciesraxperfamily = True
    self.speciesraxslowdl = True
    self.speciesraxslowdtl = True
    self.phyldog = True
    self.guenomu = False
    self.analyze = True
    self.pargenes_starting_trees = 2
    self.pargenes_bootstrap_trees = 2

  def disable_all(self):
    self.pargenes = False
    self.orthogenerax = False
    self.concatenation_min = False
    self.concatenation_max = False
    self.stag = False
    self.duptree = False
    self.fastrfs = False
    self.phyldog = False
    self.duptree = False
    self.njrax = False    
    self.njst = False    
    self.astrid = False
    self.astral = False
    self.astralpro = False
    self.speciesraxfastdl = False
    self.speciesraxfastdtl = False
    self.speciesraxprune = False
    self.speciesraxperfamily = False
    self.speciesraxslowdl = False
    self.speciesraxslowdtl = False
    self.guenomu = False
    #self.analyze = False
  
  def enable_fast_methods(self):
    self.disable_all()
    self.pargenes = True
    self.orthogenerax = True
    self.concatenation_min = True
    self.concatenation_max = True
    self.stag = True
    self.duptree = True
    self.fastrfs = True
    self.njrax = True   
    self.njst = True   
    self.astrid = True
    self.astralpro = True
    self.speciesraxfastdl = True
    self.speciesraxfastdtl = True

  def launch_reference_methods(self, datadir, subst_model, cores, launch_mode):
    misc_dir = fam.get_misc_dir(datadir)
    run_filter_file_path = os.path.join(misc_dir, str(uuid.uuid4()))
    pickle.dump(self, open(run_filter_file_path, "wb"))
    command = []
    command.append("python")
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
    if (launch_mode != "normal"):
      self.launch_reference_methods(datadir, subst_model, cores, launch_mode)
      return
    print("*************************************")
    print("Run tested species inference tools for dataset " + datadir)
    print("*************************************")
    if (self.generate):
      if (not os.path.isdir(datadir)):
        dataset = os.path.basename(datadir)
        fam_data.generate_dataset(dataset)

    if (len(datadir.split("/")) == 1):
      datadir = fam.get_datadir(datadir) 
    save_sdtout = sys.stdout
    redirected_file = os.path.join(datadir, "runs", "logs_run_all_species." + subst_model + ".txt")
    print("Redirected logs to " + redirected_file)
    sys.stdout.flush()
    sys.stdout = open(redirected_file, "w")
    if (self.pargenes):
      printFlush("Run pargenes...")
      sys.stdout.flush()
      raxml.run_pargenes_and_extract_trees(datadir, subst_model, self.pargenes_starting_trees, self.pargenes_bootstrap_trees, cores, "pargenes", True)
      sys.stdout.flush()
    if (self.stag):
      printFlush("Run Stag")
      try:
        run_stag.run_stag(datadir, subst_model)
      except Exception as exc:
        printFlush("Failed running STAG\n" + str(exc))
    if (self.duptree):
      printFlush("Run Duptree")
    try:
      run_duptree.run_duptree(datadir, subst_model)
    except Exception as exc:
      printFlush("Failed running DupTree\n" + str(exc))
    if (self.fastrfs):
      printFlush("Run FastRFS")
      try:
        run_fastrfs.run_fastrfs(datadir, subst_model)
      except Exception as exc:
        printFlush("Failed running FastRFS\n" + str(exc))
    if (self.njrax):
      printFlush("Run NJrax")
      try:
        run_njrax.run_njrax(datadir, "NJst", subst_model,cores)
      except Exception as exc:
        printFlush("Failed running NJrax\n" + str(exc))
    if (self.njst):
      printFlush("Run NJst")
      try:
        run_njst.run_njst(datadir, "raxml-ng", subst_model, "original")
        run_njst.run_njst(datadir, "raxml-ng", subst_model, "reweighted")
      except Exception as exc:
        printFlush("Failed running NJst\n" + str(exc))
    if (self.astrid):
      printFlush("Run Astral")
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
      try:
        run_astral_pro.run_astralpro(datadir, "raxml-ng", subst_model)
      except Exception as exc:
        printFlush("Failed running Astral-pro\n" + str(exc))
    if (self.concatenation_min):
      printFlush("Run concatenation-min")
      try:
        run_concatenation.run_concatenation(datadir, subst_model,False,cores)
      except Exception as exc:
        printFlush("Failed running concatenation-min\n" + str(exc))
    if (self.concatenation_max):
      printFlush("Run concatenation-max")
      try:
        run_concatenation.run_concatenation(datadir, subst_model,True,cores)
      except Exception as exc:
        printFlush("Failed running concatenation-min\n" + str(exc))
    if (self.speciesraxfastdl or self.speciesraxfastdtl):
      printFlush("Run SpeciesRaxFast")
      try:
        run_speciesrax.run_speciesrax_on_families(datadir, subst_model, cores, dl = self.speciesraxfastdl, dtl = self.speciesraxfastdtl, slow = False, strategy = "SPR")
        #run_speciesrax.run_speciesrax_on_families(datadir, subst_model, cores, dl = False, dtl = self.speciesraxfastdtl, slow = False, strategy = "TRANSFERS")
        run_speciesrax.run_speciesrax_on_families(datadir, subst_model, cores, dl = False, dtl = self.speciesraxfastdtl, slow = False, strategy = "HYBRID")
      except Exception as exc:
        printFlush("Failed running speciesrax\n" + str(exc))
    if (self.speciesraxprune):
      printFlush("Run SpeciesRaxPrune")
      try:
        dataset = os.path.basename(datadir)
        run_speciesrax.run_speciesrax_instance(dataset, "raxml-ng", True, "speciesrax-prune", subst_model, False, "SPR", cores, ["--prune-species-tree", "--per-family-rates"])
      except Exception as exc:
        printFlush("Failed running speciesrax prune\n" + str(exc))
    if (self.speciesraxperfamily):
      printFlush("Run SpeciesRaxPerFamily")
      try:
        dataset = os.path.basename(datadir)
        run_speciesrax.run_speciesrax_instance(dataset, "raxml-ng", True, "speciesrax-per-family", subst_model, False, "SPR", cores, ["--per-family-rates"])
      except Exception as exc:
        printFlush("Failed running speciesrax per family\n" + str(exc))

    if (self.speciesraxslowdl or self.speciesraxslowdtl):
      printFlush("Run SpeciesRaxSlow")
      try:
        run_speciesrax.run_speciesrax_on_families(datadir, subst_model, cores, dl = self.speciesraxslowdl, dtl = self.speciesraxslowdtl, slow = True)
      except Exception as exc:
        printFlush("Failed running speciesrax\n" + str(exc))
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
    if (self.orthogenerax):
      printFlush("Run OrthoGeneRax")
      try:
        #run_orthogenerax.run_orthogenerax(datadir, subst_model, "true", False, cores)
        run_orthogenerax.run_orthogenerax(datadir, subst_model, "njrax-NJst", False, cores)
        #run_orthogenerax.run_orthogenerax(datadir, subst_model, "speciesrax-dtl-raxml-HYBRID", False, cores)
      except Exception as exc:
        printFlush("Failed running orthogenerax\n" + str(exc))

    if (self.analyze):
      printFlush("Run analyze...")
      sys.stdout.flush()
      try:
        species_analyze.analyze(datadir)
      except Exception as exc:
        printFlush("Failed running analyze\n" + str(exc))
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


import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import fam_data
import saved_metrics
import run_all_species
import generate_families_with_subsampling
from run_all_species import SpeciesRunFilter



launch_mode = "normald"

def run_subsampled(initial_dataset, subst_model, cores, sampling_ratio, replicates, run_filter):
  initial_datadir = fam.get_datadir(initial_dataset)
  generate_families_with_subsampling.generate(initial_datadir, sampling_ratio, replicates)
  for replicate in range(0, replicates):
    dataset_dir = os.path.normpath(initial_datadir)
    dataset_dir += "_subsample" + str(sampling_ratio)
    dataset_dir += "_rep" + str(replicate)
    if (not os.path.isdir(dataset_dir)):
      print("Error: directory " + dataset_dir + " does not exist")
      sys.exit(1)
    run_filter.run_reference_methods(dataset_dir, subst_model, cores, launch_mode)



run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.duptree = True
run_filter.njrax = True
run_filter.cherry = True
run_filter.astralpro = True
#run_filter.astralpromultiple = True
#run_filter.speciesraxprune = True
#run_filter.speciesrax = True
run_filter.speciesraxperfamily = True
#run_filter.stag = True
#run_filter.astrid = True
#run_filter.orthogenerax = True
#run_filter.generaxselect = True
#run_njst = True
run_filter.cleanup = True

prot_model = "true"
dna_model = "true"

replicates = 50
if (True):
  run_subsampled("cyano_empirical", prot_model, 40, 0.01, replicates, run_filter)
  run_subsampled("ensembl_98_ncrna_primates", dna_model, 40, 0.01, replicates, run_filter)
  run_subsampled("cyano_empirical", prot_model, 40, 0.02, replicates, run_filter)
  run_subsampled("ensembl_98_ncrna_primates", dna_model, 40, 0.02, replicates, run_filter)
  run_subsampled("cyano_empirical", prot_model, 40, 0.05, replicates, run_filter)
  run_subsampled("ensembl_98_ncrna_primates", dna_model, 40, 0.05, replicates, run_filter)
  run_subsampled("cyano_empirical", prot_model, 40, 0.1, replicates, run_filter)
  run_subsampled("ensembl_98_ncrna_primates", dna_model, 40, 0.1, replicates, run_filter)




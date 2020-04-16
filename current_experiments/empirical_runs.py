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
cores = 40

run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.duptree = True
run_filter.njrax = True
run_filter.cherry = True
run_filter.astralpro = True
#run_filter.astralpromultiple = True
#run_filter.speciesraxprune = True
run_filter.speciesrax = True
run_filter.speciesraxperfamily = True
run_filter.stag = True
run_filter.orthogenerax = True
#run_filter.generaxselect = True
#run_njst = True
run_filter.cleanup = True

    
    
datasets = []
#datasets.append(("ensembl_98_ncrna_primates", "true"))
datasets.append(("ensembl_98_ncrna_mammals", "true"))
datasets.append(("ensembl_98_ncrna_sauropsids", "true"))
datasets.append(("cyano_empirical", "true"))

for dataset in datasets:
  datadir = fam.get_datadir(dataset[0])
  subst_model = dataset[1]
  run_filter.run_reference_methods(datadir, subst_model, cores, launch_mode)



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
import plot_speciesrax
import generate_families_with_grove as grove

launch_mode = "normald"
cores = 40


generate = True

run_filter = SpeciesRunFilter()
run_filter.disable_all()
run_filter.pargenes = True
run_filter.njrax = True
run_filter.cherry = True
run_filter.njst = True
run_filter.astralmulti = True
run_filter.fastmulrfs = True
run_filter.minibme = True
run_filter.minibmepruned = True
run_filter.astrid = True 
run_filter.astral = True

subst_model = "GTR+G"
minmiss = 0.0
maxmiss = 0.3
maxnodes = 200
minfam = 10
seeds = range(1000, 1100)

for seed in seeds:
  datadir = grove.get_output_dir(minmiss, maxmiss, maxnodes, minfam, seed)
  if (generate):
    grove.generate(minmiss, maxmiss, maxnodes, minfam, seed)
  print(datadir)
  assert(os.path.isdir(datadir))
  run_filter.run_reference_methods(datadir, subst_model, cores, launch_mode)







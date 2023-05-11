import subprocess
import os
import sys
import re
sys.path.insert(0, os.path.join("scripts"))
sys.path.insert(0, os.path.join("tools", "families"))
import fam
from run_all_species import SpeciesRunFilter
import generate_families_with_simphy as simphy
import experiments as exp
import shutil
import extract

def generate_ils():
  parameters =  simphy.SimphyParameters()
  parameters.prefix = "trento"
  parameters.tag = "incongruence"
  parameters.speciations_per_year = 0.000000005
  parameters.extinctions_per_year = 0.00000000
  parameters.species_taxa = 15 
  parameters.families_number = 20
  parameters.bl = 10.0    
  parameters.loss_rate = 0.0
  parameters.dup_rate = 0.0
  parameters.transfer_rate = 0.0
  parameters.gene_conversion_rate = 0.0
  parameters.sites = 1000
  parameters.model = "GTR"
  parameters.seed = 1002 # 0.089
  parameters.distance_hgt = False
  parameters.population = 100000000
  parameters.miss_species = 0.0
  parameters.miss_fam = 0.0
  parameters.subst_heterogeneity = False
  simoutput = simphy.generate_from_parameters(parameters, exp.families_datasets_root_no_fast)
  output = os.path.join(exp.families_datasets_root_no_fast, "trento_incongruence")
  families = []
  families.append("family_04")
  families.append("family_17")
  families.append("family_12")
  extract.extract(simoutput, output, families) 


def generate_tree_error():
  parameters =  simphy.SimphyParameters()
  parameters.prefix = "trento"
  parameters.tag = "error"
  parameters.speciations_per_year = 0.000000005
  parameters.extinctions_per_year = 0.00000000
  parameters.species_taxa = 15 
  parameters.families_number = 50
  parameters.bl = 0.5
  parameters.loss_rate = 0.0
  parameters.dup_rate = 0.0
  parameters.transfer_rate = 0.0
  parameters.gene_conversion_rate = 0.0
  parameters.sites = 50
  parameters.model = "GTR"
  parameters.seed = 1000
  parameters.distance_hgt = False
  parameters.population = 10
  parameters.miss_species = 0.0
  parameters.miss_fam = 0.0
  parameters.subst_heterogeneity = False
  simoutput = simphy.generate_from_parameters(parameters, exp.families_datasets_root_no_fast)
  output = os.path.join(exp.families_datasets_root_no_fast, "trento_error")
  families = []
  for family in fam.get_families_list(simoutput):
    families.append(family)
  extract.extract(simoutput, output, families) 

if (__name__ == "__main__"): 
  generate_ils()
  generate_tree_error()


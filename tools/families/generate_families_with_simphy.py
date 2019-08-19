import os
import sys
import subprocess
import shutil
import fam
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "phyldog"))
sys.path.insert(0, os.path.join("tools", "trees"))
import experiments as exp

def build_config_file(output_dir):
  species_taxa = 10
  families_number = 5
  config_file = os.path.join(output_dir, "simphy_config.txt")
  with open(config_file, "w") as writer:
    writer.write("// SPECIES TREE\n")
    writer.write("-RS 1 // number of replicates\n")
    writer.write("-sb f:0.00001 // speciations per year\n")
    writer.write("-sl f:" + str(species_taxa) + " // species taxa\n")

    writer.write("// POPULATION\n")
    writer.write("-SP f:1000 // population size\n")

    writer.write("// LOCUS\n")
    writer.write("-rl f:" + str(families_number) + " // locus (gene family) per replicate\n")

    writer.write("// GENERAL\n")
    writer.write("-cs 42 // random seed\n")
    writer.write("-O " + str(output_dir) + " // output directory\n")
    writer.write("-OM 1 // output the mappings\n")
    writer.write("-OC 1 // log the configuration file\n")

  return config_file

def run_simphy(config_file):
  commands = []
  commands.append(exp.simphy_exec)
  commands.append("-I")
  commands.append(config_file)
  subprocess.check_call(commands)

def export_to_family(output_dir):
  pass

output_dir = "simphy_test"
exp.reset_dir(output_dir)
config_file = build_config_file(output_dir)
run_simphy(config_file)
export_to_family(output_dir)







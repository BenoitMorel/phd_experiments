import os
import sys
import subprocess
sys.path.insert(0, 'scripts')
import experiments as exp

def generate_zombi_species(species, output):
  parameters_dir = os.path.join(output, "parameters")
  species_parameters_file = os.path.join(output, "SpeciesTreeParameters.tsv")
  with open(species_parameters_file, "w") as writer:
    writer.write("SPECIATION f:0.1\n")
    writer.write("EXTINCTION f:0.02\n")
    writer.write("STOPPING_RULE 1\n")
    writer.write("TOTAL_LINEAGES " + str(species) + "\n")
    writer.write("TOTAL_TIME 1\n")
    writer.write("MIN_LINEAGES 1\n")
    writer.write("MAX_LINEAGES 10000\n")
    writer.write("VERBOSE 1\n")
    writer.write("TURNOVER F:0.0002\n")
    writer.write("LINEAGE_PROFILE 100-100;300-15000;500-50\n")
  command = []
  command.append("python3")
  command.append(exp.zombi_script)
  command.append("T")
  command.append(species_parameters_file)
  command.append(output)
  subprocess.check_call(command)

def generate_zombi_genome(families, dupRate, lossRate, transferRate, output):
  parameters_dir = os.path.join(output, "parameters")
  genome_parameters_file = os.path.join(output, "GenomeTreeParameters.tsv")
  with open(genome_parameters_file, "w") as writer:
    writer.write("DUPLICATION f:" + str(dupRate) + "\n")
    writer.write("TRANSFER f:" + str(transferRate) + "\n")
    writer.write("LOSS f:" + str(lossRate) + "\n")
    writer.write("INVERSION f:0\n")
    writer.write("TRANSLOCATION f:0\n")
    writer.write("ORIGINATION f:0\n")
    writer.write("DUPLICATION_EXTENSION f:1\n")
    writer.write("TRANSFER_EXTENSION f:1\n")
    writer.write("LOSS_EXTENSION f:1\n")
    writer.write("INVERSION_EXTENSION f:0.05\n")
    writer.write("TRANSLOCATION_EXTENSION f:0.3\n")
    writer.write("REPLACEMENT_TRANSFER 1\n")
    writer.write("INITIAL_GENOME_SIZE " + str(families) + "\n")
    writer.write("MIN_GENOME_SIZE 10\n")
    writer.write("GENE_LENGTH f:100\n")
    writer.write("INTERGENE_LENGTH 100\n")
    writer.write("PSEUDOGENIZATION 0.5\n")
    writer.write("####OUTPUT\n")
    writer.write("EVENTS_PER_BRANCH 1\n")
    writer.write("PROFILES 1\n")
    writer.write("GENE_TREES 1\n")
    writer.write("RECONCILED_TREES 0\n")
    writer.write("VERBOSE 1\n")

  command = []
  command.append("python3")
  command.append(exp.zombi_script)
  command.append("G")
  command.append(genome_parameters_file)
  command.append(output)
  subprocess.check_call(command)
  
def generate_zombi_sequence(sites, output):
  parameters_dir = os.path.join(output, "parameters")
  sequence_parameters_file = os.path.join(output, "SequenceTreeParameters.tsv")
  states = ["A", "C", "G", "T"]
  with open(sequence_parameters_file, "w") as writer:
    writer.write("SCALING 1\n")
    writer.write("SEQUENCE_SIZE " + str(sites) + "\n")
    writer.write("SEQUENCE codon\n")
    for s1 in states:
      for s2 in states:
        if (s1 != s2):
          writer.write(s1 +  s2 + " 1.0\n")
    for s in states:
      writer.write(s + " 0.25\n")
    writer.write("CODON_MODEL MG\n")
    writer.write("ALPHA 1.0\n")
    writer.write("BETA 0.5\n")
    writer.write("KAPPA 1.0\n")
    writer.write("VERBOSE 1\n")
  
  command = []
  command.append("python3")
  command.append(exp.zombi_script)
  command.append("S")
  command.append(sequence_parameters_file)
  command.append(output)
  subprocess.check_call(command)


def generate_zombi(species, families, sites, dupRate, lossRate, transferRate, output):
  os.makedirs(output)
  parameters_dir = os.path.join(output, "parameters")
  os.makedirs(parameters_dir)
  generate_zombi_species(species, output) 
  generate_zombi_genome(families, dupRate, lossRate, transferRate, output) 
  generate_zombi_sequence(sites, output)

if (len(sys.argv) != 8):
  print("Syntax: python generate_zombi.py species families sites dupRate lossRate transferRate output")
  sys.exit(1)

species = int(sys.argv[1])
families = int(sys.argv[2])
sites = int(sys.argv[3])
dupRate = float(sys.argv[4])
lossRate = float(sys.argv[5])
transferRate = float(sys.argv[6])
output = sys.argv[7]
generate_zombi(species, families, sites, dupRate, lossRate, transferRate, output)

import sys
import os
import math

def get_model_free_parameters(model):
  free_parameters = 0
  if ("+G" in model):
    free_parameters += 1
  if ("+I" in model):
    free_parameters += 1
  if ("+F" in model):
    free_parameters += 19
  if (model.startswith("LG4M")):
    free_parameters += 1
  if (model.startswith("LG4X")):
    free_parameters += 6
  #print(model + ": " + str(free_parameters))
  return free_parameters

if (len(sys.argv) != 2):
  print("usage: python raxml_free_parameters.py raxml_logs")
  sys.exit(1)
input_part = sys.argv[1]

models_free_parameters = 0
free_parameters = 0
brlen = ""
sites = 0
taxa = 0
branches_count = 0
partitions_count = 0
log_likelihood = 0.0

lines = open(input_part).readlines()
dico_partitions = {}
for idx, line in enumerate(lines):
  if (line.startswith("Model:")):
# first check that we did not treat this partition
# already (in case of a restart, raxml logs are doubled)
    line_partition = lines[idx-1]
    if (not line_partition.startswith("Partition")):
      print("Error when reading partition")
      sys.exit(1)
    partition = line_partition.split(" ")[2]
    if (partition in dico_partitions):
      continue
    dico_partitions[partition] = partition
    model = line.split(" ")[1][:-1]
    models_free_parameters += get_model_free_parameters(model)
    partitions_count += 1
  elif ("branch lengths:" in line):
    brlen = line.split(" ")[4]
  elif ("Loaded align" in line):
    sites = int(line.split(" ")[7])
    taxa = int(line.split(" ")[4])
    branches_count = 2 * taxa - 2
  elif (line.startswith("Final LogLikelihood:")):
    log_likelihood = float(line.split(" ")[2][:-1])
  elif ("Tree #1, final logLikelihood:" in line): # --eval case
    log_likelihood = float(line.split(" ")[5][:-1])


free_parameters = models_free_parameters
if (brlen == "linked"):
  free_parameters += branches_count
elif (brlen == "unlinked"):
  free_parameters += branches_count * partitions_count
elif (brlen == "proportional"):
  free_parameters += branches_count + partitions_count - 1
else:
  print("Unrecognized brlen " + brlen + ". Aborting.")
  sys.exit(1)

sample_size = sites
if (free_parameters >= sample_size):
  print("more free parameters than samples! I will not compute the Akaike scores!")
  sys.exit(0)

AIC = 2.0 * free_parameters - 2.0 * log_likelihood
BIC = math.log(sample_size) * free_parameters - 2.0 * log_likelihood
AICc = AIC + 2.0 * free_parameters * (free_parameters + 1) / (sample_size - free_parameters - 1) 

print("")
print("Brlen mode: " + brlen)
print("Partitions: " + str(partitions_count))
print("Sites: " + str(sites))
print("Taxa: " + str(taxa))
print("Branches: " + str(branches_count))
print("Models free parameters: " +  str(models_free_parameters))
print("Total free parameters: " + str(free_parameters))
print("Log likelihood: " + str(log_likelihood))
print("")
print("AIC:  " + str(AIC))
print("AICc: " + str(AICc))
print("BIC:  " + str(BIC))



import sys
import os
import math

def get_model_free_rates(model):
  free_rates = 0
  if ("+G" in model):
    free_rates += 1
  if ("+I" in model):
    free_rates += 1
  if ("+F" in model):
    free_rates += 19
  if (model.startswith("LG4M")):
    free_rates += 1
  if (model.startswith("LG4X")):
    free_rates += 6
  #print(model + ": " + str(free_rates))
  return free_rates

if (len(sys.argv) != 2):
  print("usage: python raxml_free_rates.py raxml_logs")
  sys.exit(1)
input_part = sys.argv[1]

models_free_rates = 0
free_rates = 0
brlen = ""
sites = 0
taxa = 0
branches_count = 0
partitions_count = 0
log_likelihood = 0.0

lines = open(input_part).readlines()
for line in lines:
  if (line.startswith("Model:")):
    model = line.split(" ")[1][:-1]
    models_free_rates += get_model_free_rates(model)
    partitions_count += 1
  elif ("branch lengths:" in line):
    brlen = line.split(" ")[4]
  elif ("Loaded align" in line):
    sites = int(line.split(" ")[7])
    taxa = int(line.split(" ")[4])
    branches_count = 2 * taxa - 2
  elif (line.startswith("Final LogLikelihood:")):
    log_likelihood = float(line.split(" ")[2][:-1])


free_rates = models_free_rates
if (brlen == "linked"):
  free_rates += branches_count
elif (brlen == "unlinked"):
  free_rates += branches_count * partitions_count
elif (brlen == "scaled"):
  free_rates += branches_count + partitions_count - 1
else:
  print("Unrecognized brlen " + brlen + ". Aborting.")
  sys.exit(1)

AIC = 2.0 * free_rates - 2.0 * log_likelihood
BIC = math.log(sites) * free_rates - 2.0 * log_likelihood
AICc = AIC + 2.0 * free_rates * (free_rates + 1) / (sites - free_rates - 1) 

print("")
print("Brlen mode: " + brlen)
print("Partitions: " + str(partitions_count))
print("Sites: " + str(sites))
print("Taxa: " + str(taxa))
print("Branches: " + str(branches_count))
print("Models free rates: " +  str(models_free_rates))
print("Total free rates: " + str(free_rates))
print("Log likelihood: " + str(log_likelihood))
print("")
print("AIC:  " + str(AIC))
print("AICc: " + str(AICc))
print("BIC:  " + str(BIC))



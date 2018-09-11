import sys
import os


if (len(sys.argv) != 2):
  print("usage: python compute_free_rates.py input_part")
  sys.exit(1)
input_part = sys.argv[1]


total_free_rates = 0
lines = open(input_part).readlines()
for line in lines:
  free_rates = 0
  model = line.split(",")[0]
  if ("+G" in model):
    free_rates += 1
  if ("+I" in model):
    free_rates += 1
  if ("+F" in model):
    free_rates += 1
  if (model.startswith("LG4M")):
    free_rates += 1
  if (model.startswith("LG4X")):
    free_rates += 6
  print(model + " " + str(free_rates))
  total_free_rates += free_rates
print("Total number of free rates: " + str(total_free_rates))


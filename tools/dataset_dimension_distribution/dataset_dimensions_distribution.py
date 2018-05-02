import os
import sys
import seaborn
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


if (len(sys.argv) != 3):
  print("Syntax error. Usage: python script.py results_directory plot_name")

output_dir = sys.argv[1]
plot_name = sys.argv[2]

parsing_results = os.path.join(output_dir, "parse_run", "results")
msa_names = os.listdir(parsing_results)

print("msas: " + str(msa_names))
sites_array = []
taxa_array = []
for name in msa_names:
  log_file = os.path.join(parsing_results, name, name + ".raxml.log")
  unique_sites = 0
  taxa = 0
  try:
    lines = open(log_file).readlines()
  except:
    continue
  for line in lines:
    if "Alignment comprises" in line:
      unique_sites = int(line.split(" ")[5])
    if "taxa" in line:
      taxa = int(line.split(" ")[4])
  if (unique_sites * taxa == 0):
    continue
  sites_array.append(unique_sites)
  taxa_array.append(taxa)
  print(str(taxa) + " " + str(unique_sites))

default = range(0, len(taxa_array))

fig, ax = plt.subplots()
ax.scatter(taxa_array, sites_array, marker='x')
ax.grid()
plt.xlabel("taxa")
plt.ylabel("sites")
fig.savefig(plot_name + "_sites_taxa.png")

plt.clf()
fig, ax = plt.subplots()
ax.scatter(default, sites_array, marker='x')
plt.ylabel("sites")
ax.grid()
fig.savefig(plot_name + "_sites.png")

plt.clf()
fig, ax = plt.subplots()
ax.scatter(default, taxa_array, marker='x')
plt.ylabel("taxa")
ax.grid()
fig.savefig(plot_name + "_taxa.png")



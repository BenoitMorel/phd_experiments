import os
import sys
import seaborn
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

if (len(sys.argv) != 3 or len(sys.argv) != 5):

  print("Syntax error. Usage: python script.py results_directory plot_name [maxX maxY]")

output_dir = sys.argv[1]
plot_name = sys.argv[2]

maxX = 0
maxY = 0
if (len(sys.argv) == 5):
  maxX = int(sys.argv[3])
  maxY = int(sys.argv[4])

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
ax.scatter(sites_array, taxa_array, marker='x')
ax.grid()
plt.xlabel("sites")
plt.ylabel("taxa")
axes = plt.gca()
if (maxX != 0):
  axes.set_xlim([0,maxX])
if (maxY != 0):
  axes.set_ylim([0,maxY])
fig.savefig(plot_name + "_sites_taxa.png")

plt.clf()
fig, ax = plt.subplots()
ax.scatter(default, sites_array, marker='x')
plt.ylabel("sites")
ax.grid()
if (maxX != 0):
  axes.set_xlim([0,maxX])
if (maxY != 0):
  axes.set_ylim([0,maxY])
fig.savefig(plot_name + "_sites.png")

plt.clf()
fig, ax = plt.subplots()
ax.scatter(default, taxa_array, marker='x')
plt.ylabel("taxa")
ax.grid()
if (maxX != 0):
  axes.set_xlim([0,maxX])
if (maxY != 0):
  axes.set_ylim([0,maxY])
fig.savefig(plot_name + "_taxa.png")

print("result in " + plot_name + "_sites_taxa")

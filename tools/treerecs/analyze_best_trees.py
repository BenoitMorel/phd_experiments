import os
import sys

import matplotlib.pyplot as plt
#import seaborn as sns



class TreeEntry(object):
  def __init__(self):
    self.family = ""
    self.tree = ""
    self.threshold = 0
    self.ale_ll = 0.0
    self.pll_ll = 0.0
    self.joint_ll = 0.0

  def __str__(self):
    return "(" + str(self.family) + "," + str(self.threshold) + "," + str(self.ale_ll)+ "," + str(self.pll_ll) + "," + str(self.joint_ll)


class AllTreeEntries(object):
  def __init__ (self, treerecs_output):
    self.all_entries = {}
    self.best_entries = {}
    self.per_threshold_trees = {}
    self.__fill(treerecs_output)

  def get_best_thresholds_histogram(self):
    histo = {}
    for family in self.best_entries:
      threshold = self.best_entries[family].threshold
      if (not threshold in histo):
        histo[threshold] = 0
      histo[threshold] += 1
    return histo

  def __add_entry(self, entry):
    family = entry.family
    threshold = entry.threshold
    family_entries = None
    if (not threshold in self.per_threshold_trees):
      self.per_threshold_trees[threshold] = []
    self.per_threshold_trees[threshold].append(entry.tree)
    if (not family in self.all_entries):
      family_entries = {}
      self.all_entries[family] = family_entries
      self.best_entries[family] = entry
    else:
      family_entries = self.all_entries[family]
      if (self.best_entries[family].joint_ll <= entry.joint_ll):
        self.best_entries[family] = entry
    family_entries[threshold] = entry
    #print("Entry: " + str(entry))
    #print("Best : " + str(self.best_entries[family]))

  def __add_entry_from_line(self, line, tree):
    split = line.split(" ")
    entry = TreeEntry()
    entry.tree = tree
    for index, val in enumerate(split):
        if (val == "family"):
          entry.family = split[index + 1]
        elif (val == "contraction" and split[index + 1] == "threshold"):
          threshold = split[index + 3][:-1]
          if (threshold == "no"):
              threshold = "0.0"
          entry.threshold = float(threshold)
        elif (val == "logLk"):
          if (split[index - 1] == "ALE"):
            entry.ale_ll = float(split[index + 2][:-1])
          elif (split[index - 1] == "libpll"):
            entry.pll_ll = float(split[index + 2][:-1])
    entry.joint_ll = entry.pll_ll + entry.ale_ll
    self.__add_entry(entry)

  def __fill(self, treerecs_output):
    lines = open(treerecs_output).readlines()
    previous_line = None
    for line in lines:
      if (not line.startswith(">")):
        self.__add_entry_from_line(previous_line, line[:-1])
      previous_line = line

  def save_per_threshold_trees(self, path, basename):
    for threshold in self.per_threshold_trees:
      trees = self.per_threshold_trees[threshold]
      name = basename + str(int(threshold * 100)) + ".txt"
      output = os.path.join(path, name)
      with open(output, "w") as f:
        for tree in trees:
          f.write(tree + "\n")

def plot_histogram(histogram, output_file, xlabel, ylabel):
  fig, ax = plt.subplots()
  x = histogram.keys()
  y = histogram.values()
  plt.bar(x, y, 0.04)
  for i,j in zip(x,y):
        ax.annotate(str(j),xy=(i,j))
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  fig.savefig(output_file)



if (len(sys.argv) != 3):
  print("Syntax: python analyze_best_trees.py treerecs_output analyze_output_dir")
  sys.exit(1)

treerecs_output = sys.argv[1]
output = sys.argv[2]
trees_output = os.path.join(output, "trees")
all_entries = AllTreeEntries(treerecs_output)

try:
  os.makedirs(trees_output)
except:
  pass
histo = all_entries.get_best_thresholds_histogram()
print(sorted(histo.items(), key=lambda x: x[0]))
plot_histogram(histo, os.path.join(output, "threshold_histogram.png"), "threshold", "best trees number")
all_entries.save_per_threshold_trees(trees_output, "trees_for_threshold_")



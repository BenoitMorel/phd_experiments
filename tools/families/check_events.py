import sys
import os
from ete3 import Tree
import numpy
import math
sys.path.insert(0, os.path.join("tools", "trees"))
from read_tree import read_tree
from rf_distance import ete3_rf
from rf_distance import get_relative_rf
from rf_distance import get_rf

def get_nodes(tree1):
  rf = ete3_rf(tree1, tree1)
  return rf[1]

def hamming(v1, v2):
  dist = 0
  zipVector = zip(v1, v2)
  for member in zipVector:
    if (member[0] != member[1]):
      dist += 1
  return float(dist) / float(len(v1))

def squared(v1, v2):
  quar_distance = 0
  zipVector = zip(v1, v2)
  for member in zipVector:
    quar_distance += (member[1] - member[0]) ** 2
  return float(quar_distance) / float(len(v1) ** 2)

def euclidean (v1, v2):
  quar_distance = 0
  zipVector = zip(v1, v2)
  for member in zipVector:
    quar_distance += (member[1] - member[0]) ** 2
  return math.sqrt(quar_distance) / float(len(v1))


def fstr(a):
  return "{0:.4f}".format(a)

def align(method):
  return method + " " * (15 - len(method))

def printDistances(events, method, event_type):
  v1 = events[method][event_type]
  v2 = events["True"][event_type]
  if (len(v1) != len(v2) or len(v1) == 0):
    return 
  print("- " + align(method) + ":  euclidean = " + fstr(euclidean(v1, v2)) + "  hamming = " + fstr(hamming(v1, v2)))

def analyze_events(dataset_dir, analyze_dir):
  event_types = ["S", "D", "T"]
  methods = ["True", "Treerecs", "Phyldog", "JointSearch"]
  suffixes = {}
  prefixes = {}
  prefixes["True"] = dataset_dir
  suffixes["True"] = "trueEvents.txt"
  prefixes["Treerecs"] = dataset_dir
  suffixes["Treerecs"] =  "treerecsEvents.txt"
  prefixes["Phyldog"] = dataset_dir
  suffixes["Phyldog"] = "phyldogEvents.txt" 
  prefixes["JointSearch"] = os.path.join(analyze_dir, "results")
  suffixes["JointSearch"] = "jointsearch.events"

  events = {}
  for method in methods:
    events[method] = {}
    for event_type in event_types:
      events[method][event_type] = [] 

  for msa in os.listdir(dataset_dir): 
    for method in methods:
      event_file = None
      try:
        events_file = os.path.join(prefixes[method], msa, suffixes[method])
        events_lines = open(events_file).readlines()
        if (len(events_lines) == 0 or len(events_lines[0]) == 0):
          continue
      except:
        continue
      events[method]["S"].append(int(events_lines[0].split(":")[1][:-1]))
      if (method == "JointSearch"):
        events[method]["D"].append(int(events_lines[2].split(":")[1][:-1]))
        events[method]["T"].append(int(events_lines[3].split(":")[1][:-1]) + int(events_lines[4].split(":")[1][:-1]))
      else:
        events[method]["D"].append(int(events_lines[1].split(":")[1][:-1]))
        events[method]["T"].append(int(events_lines[2].split(":")[1][:-1]))

  transferPresent = (sum(events["True"]["T"]) != 0) or (sum(events["JointSearch"]["T"]) != 0)
 
  if (False):
    print("Duplications: ")
    for method in methods:
      if (len(events[method]["D"]) == 0):
        continue
      print(method + " " + str(events[method]["D"]))
    print("")

    if (transferPresent):
      print("Transfers: ")
      for method in methods:
        if (len(events[method]["T"]) == 0):
          continue
        print(method + " " + str(events[method]["T"]))
    print("")

  print("Duplication event count vectors (normalized distances with true vectors)")
  for method in methods:
    if (method == "True"):
      continue
    printDistances(events, method, "D")
  print("")

  if (transferPresent):
    print("Transfer event count vectors (normalized distances with true vectors)")
    for method in methods:
      if (method == "True"):
        continue
      printDistances(events, method, "T")
  print("")
  

if __name__ == '__main__':
  if (len(sys.argv) < 3):
    print("Syntax: families_dir analyze_dir")
    exit(1)
  print(" ".join(sys.argv))
  dataset_dir = sys.argv[1]
  analyze_dir = sys.argv[2]
  analyze_events(dataset_dir, analyze_dir)




import sys
import os
from ete3 import Tree
import numpy
import math

def read_tree(tree):
  lines = open(tree).readlines()
  for line in lines:
    if (not line.startswith(">")):
      return Tree(line, format=1)
  return None

def get_rf(tree1, tree2):
  return tree1.robinson_foulds(tree2, unrooted_trees=True)[0]

def get_relative_rf(tree1, tree2):
  rf = tree1.robinson_foulds(tree2, unrooted_trees=True)
  return float(rf[0]) / float(rf[1])

def get_nodes(tree1):
  rf = tree1.robinson_foulds(tree1, unrooted_trees=True)
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

def analyse_events(dataset_dir, jointsearch_scheduler_dir):
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
  prefixes["JointSearch"] = os.path.join(jointsearch_scheduler_dir, "results")
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

def analyse(dataset_dir, jointsearch_scheduler_dir):
  analysed_msas = 0
  total_nodes_number = 0
  methods = ["Raxml-ng", "Treerecs", "Phyldog", "JointSearch"]
  methods_tree_files = {}
  methods_tree_files["Raxml-ng"] = "raxmlGeneTree.newick"
  methods_tree_files["Treerecs"] = "treerecsGeneTree.newick"
  methods_tree_files["Phyldog"] = "phyldogGeneTree.newick"
  methods_tree_files["JointSearch"] = "jointsearch.newick"
  methods_trees_number = {}
  js_dup = []
  js_loss = []
  js_trans = []
  js_initialll = []
  js_initialllrec = []
  js_initiallllibpll = []
  js_ll = []
  js_llrec = []
  js_lllibpll = []
  for m in methods:
    methods_trees_number[m] = 0
  total_rrf = {}
  total_rf = {}
  true_matches = {}
  best_tree = {}
  tca_scores = {}
  good_tca_scores = []
  bad_tca_scores = []
  for method in methods:
    total_rrf[method] = 0.0
    total_rf[method] = 0.0
    true_matches[method] = 0
    best_tree[method] = 0
  for msa in os.listdir(dataset_dir):   
    trees = {}
    family_path = os.path.join(dataset_dir, msa)
    tca_path = os.path.join(family_path, "tca.txt")
    try:
      tca_scores[msa] = float(open(tca_path).readlines()[0])
    except:
      pass
    try:
      true_tree = Tree(os.path.join(family_path, "trueGeneTree.newick"), format=1) 
    except:
      continue
    jointsearch_prefix = os.path.join(jointsearch_scheduler_dir, "results", msa)
    for method in methods:
      if (method == "JointSearch"):
        prefix = jointsearch_prefix
      else:
        prefix = family_path
      try:
        trees[method] = read_tree(os.path.join(prefix, methods_tree_files[method]))
        methods_trees_number[method] += 1
      except:
        try:
          trees[method] = Tree(os.path.join(family_path, "trueGeneTree.newick"), format=1)
        except:
          trees[method] = None
    best_rrf = 1
    rrf = {}
    rf = {}
    analysed_msas += 1
    for method in methods:
      if (trees[method] == None):
        continue
      rf_cell = trees[method].robinson_foulds(true_tree, unrooted_trees=True)
      rrf[method] = float(rf_cell[0]) / float(rf_cell[1])
      rf[method] = float(rf_cell[0])
      best_rrf = min(best_rrf, rrf[method])
      total_rrf[method] += rrf[method]
      total_rf[method] += float(rf_cell[0])
      if (rf_cell[0] == 0):
        true_matches[method] += 1
      if (method == methods[0]): #do it only once
        total_nodes_number += rf_cell[1]
    for method in methods:
      if (best_rrf == rrf[method]):
        best_tree[method] += 1
    if (rf["Raxml-ng"] == 0.0):
      good_tca_scores.append(tca_scores[msa])
    else:
      bad_tca_scores.append(tca_scores[msa])
    stats_file = os.path.join(jointsearch_prefix, "jointsearch.stats")
    with open(stats_file) as stats_reader:
      lines = stats_reader.readlines()
      js_initialll.append(float(lines[0].split(" ")[1][:-1]))
      js_initialllrec.append(float(lines[1].split(" ")[1][:-1]))
      js_initiallllibpll.append(float(lines[2].split(" ")[1][:-1]))
      js_ll.append(float(lines[3].split(" ")[1][:-1]))
      js_llrec.append(float(lines[4].split(" ")[1][:-1]))
      js_lllibpll.append(float(lines[5].split(" ")[1][:-1]))
      js_dup.append(float(lines[6].split(" ")[1][:-1]))
      js_loss.append(float(lines[7].split(" ")[1][:-1]))
      js_trans.append(float(lines[8].split(" ")[1][:-1]))
  if (analysed_msas == 0):
    print("did not manage to analyse any MSA")
    exit(1)
  
  if (len(good_tca_scores) != 0):
    print("average and min tca for raxml trees that match true trees: " + str(numpy.mean(good_tca_scores)) + " " + str(min(good_tca_scores)) + " (on " + str(len(good_tca_scores)) + " elements)")
  if (len(bad_tca_scores) != 0):

    print("average and max tca for raxml trees that DO NOT match true trees: " + str(numpy.mean(bad_tca_scores)) + " " + str(max(bad_tca_scores)) + " (on " + str(len(bad_tca_scores)) + " elements)")

  print("Rates arrays")
  print("D:")
  print(js_dup)
  print("L:")
  print(js_loss)
  print("T:")
  print(js_trans)
  print("")

  print("Number of gene families: " + str(analysed_msas))
  print("")

  print("Total initial joint likelihood: " + str(sum(js_initialll)))
  print("Total initial libpll  likelihood: " + str(sum(js_initiallllibpll)))
  print("Total initial reconciliation likelihood: " + str(sum(js_initialllrec)))
  print("")
  print("Total joint likelihood: " + str(sum(js_ll)))
  print("Total libpll  likelihood: " + str(sum(js_lllibpll)))
  print("Total reconciliation likelihood: " + str(sum(js_llrec)))
  print("")
  print("Average D=" + str(numpy.mean(js_dup)))
  print("Average L=" + str(numpy.mean(js_loss)))
  print("Average T=" + str(numpy.mean(js_trans)))
  print("")
  print("Standard deviation D=" + str(numpy.std(js_dup)))
  print("Standard deviation L=" + str(numpy.std(js_loss)))
  print("Standard deviation T=" + str(numpy.std(js_trans)))
  print("")

  print("Average (over the gene families) relative RF distance to the true trees:")
  for method in methods:
    print("- " + method + ":\t" + str(total_rrf[method] / float(analysed_msas)))
  print("")
  
  print("Normalized average (over the gene families) RF distance to the true trees:")
  for method in methods:
    print("- " + method + ":\t" + str(total_rf[method] / float(total_nodes_number)))
  print("")
  
  print("Number of gene families for which a method reaches the smallest relative RF to the true trees compared with the other methods:")
  for method in methods:
    print("- " + method + ":\t" + str(best_tree[method]) + "/" + str(analysed_msas))
  print("")
  
  print("Number of gene families for which a method finds the true tree:")
  for method in methods:
    print("- " + method + ":\t" + str(true_matches[method]) + "/" + str(analysed_msas))
  print("")

  print("Analysed trees:")
  for method in methods:
    print("- " + method + ":\t" + str(methods_trees_number[method]))
  print("")

  analyse_events(dataset_dir, jointsearch_scheduler_dir)

if __name__ == '__main__':
  if (len(sys.argv) != 3):
    print("Syntax: families_dir jointsearch_scheduler_dir")
    exit(1)
  print(" ".join(sys.argv))
  dataset_dir = sys.argv[1]
  jointsearch_scheduler_dir = sys.argv[2]
  analyse(dataset_dir, jointsearch_scheduler_dir)




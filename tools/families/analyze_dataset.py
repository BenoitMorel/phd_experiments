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

def analyse_events(dataset_dir, analyze_dir):
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
  
methods_tree_files = {}
methods_tree_files["True"] = "trueGeneTree.newick"
methods_tree_files["RAxML-NG"] = "raxmlGeneTree.newick"
methods_tree_files["Treerecs"] = "treerecsGeneTree.newick"
methods_tree_files["Phyldog"] = "phyldogGeneTree.newick"
methods_tree_files["Notung"] = "notungGeneTree.newick"
methods_tree_files["JointSearch"] = "jointsearch.newick"
methods_tree_files["GeneRax"] = "geneTree.newick"

def get_gene_tree(method, dataset_dir, analyze_dir, msa):
  if (method == "JointSearch"):
    prefix = os.path.join(analyze_dir, "results", msa)
  elif (method == "GeneRax"):
    prefix = os.path.join(analyze_dir, msa)
  else:
    prefix = os.path.join(dataset_dir, msa)
  return read_tree(os.path.join(prefix, methods_tree_files[method]))

def analyse(dataset_dir, analyze_dir, benched_method = "JointSearch"):
  analysed_msas = 0
  total_nodes_number = 0
  methods = ["True", "RAxML-NG", "Treerecs", "Phyldog", "Notung", benched_method]
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
  methods_to_compare = []
  methods_to_compare.append(("True", "RAxML-NG"))
  methods_to_compare.append(("True", "Treerecs"))
  methods_to_compare.append(("True", "Phyldog"))
  methods_to_compare.append(("True", "Notung"))
  methods_to_compare.append(("True", benched_method))
  methods_to_compare.append((benched_method, "RAxML-NG"))
  methods_to_compare.append((benched_method, "Treerecs"))
  for m in methods:
    methods_trees_number[m] = 0
  total_rrf = {}
  best_tree = {}
  for method in methods:
    best_tree[method] = 0
  for method_pair in methods_to_compare:
    method1 = method_pair[0]
    method2 = method_pair[1]
    methods_key = method1 + " - " + method2
    total_rrf[methods_key] = 0.0
  for msa in os.listdir(dataset_dir):   
    trees = {}
    family_path = os.path.join(dataset_dir, msa)
    jointsearch_prefix = os.path.join(analyze_dir, "results", msa)
    invalid_methods = []
    for method in methods:
      try: 
        trees[method] = get_gene_tree(method, dataset_dir, analyze_dir, msa)
      except:
        print("Cannot read " + method + " tree for " + msa)
        invalid_methods.append(method)
    for method in invalid_methods:
      methods.remove(method)
      methods_to_compare = [p for p in methods_to_compare if (p[0] != method and p[1] != method)]
      print("Missing tree for " + method)
      print("This method will be excluded")
    best_rrf = 1
    rrf = {}
    analysed_msas += 1
    for method_pair in methods_to_compare:
      method1 = method_pair[0]
      method2 = method_pair[1]
      methods_key = method1 + " - " + method2
      rf_cell = ete3_rf(trees[method1], trees[method2])
      if (rf_cell[1] == 0):
        print("null cell for " + methods_key + " " + msa)
        exit(1)
      rrf[methods_key] = float(rf_cell[0]) / float(rf_cell[1])
      if (method1 == "True" or method2 == "True"):
        best_rrf = min(best_rrf, rrf[methods_key])
      total_rrf[methods_key] += rrf[methods_key]
    for method_pair in methods_to_compare:
      method1 = method_pair[0]
      method2 = method_pair[1]
      methods_key = method1 + " - " + method2
      if (method1 == "True"):
        if (best_rrf == rrf[methods_key]):
          best_tree[method2] += 1
      elif (method2 == "True"):
        if (best_rrf == rrf[methods_key]):
          best_tree[method1] += 1
    stats_file = os.path.join(jointsearch_prefix, "jointsearch.stats")
    try:
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
    except:
      pass
  if (analysed_msas == 0):
    print("did not manage to analyse any MSA")
    exit(1)
  

  #print("Rates arrays")
  #print("D:")
  #print(js_dup)
  #print("L:")
  #print(js_loss)
  #print("T:")
  #print(js_trans)
  #print("")

  print("Number of gene families: " + str(analysed_msas))
  print("")

  if (len(js_ll) > 0):
    #print("Total initial joint likelihood: " + str(sum(js_initialll)))
    #print("Total initial libpll  likelihood: " + str(sum(js_initiallllibpll)))
    #print("Total initial reconciliation likelihood: " + str(sum(js_initialllrec)))
    #print("")
    print("Total joint likelihood: " + str(sum(js_ll)))
    #print("Total libpll  likelihood: " + str(sum(js_lllibpll)))
    #print("Total reconciliation likelihood: " + str(sum(js_llrec)))
    #print("")
    print("Average D=" + str(numpy.mean(js_dup)))
    print("Average L=" + str(numpy.mean(js_loss)))
    print("Average T=" + str(numpy.mean(js_trans)))
    print("")
    #print("Standard deviation D=" + str(numpy.std(js_dup)))
    #print("Standard deviation L=" + str(numpy.std(js_loss)))
    #print("Standard deviation T=" + str(numpy.std(js_trans)))
    #print("")

  print("Average (over the gene families) relative RF distances:")
  for method_pair in methods_to_compare:
    method1 = method_pair[0]
    method2 = method_pair[1]
    methods_key = method1 + " - " + method2
    print("- " + methods_key + ":\t" + str(total_rrf[methods_key] / float(analysed_msas)))
  print("")
  
  
  print("Number of gene families for which a method reaches the smallest relative RF to the true trees compared with the other methods:")
  for method in methods:
    if (method == "True"):
      continue
    print("- " + method + ":\t" + str(best_tree[method]) + "/" + str(analysed_msas))
  print("")
  

  analyse_events(dataset_dir, analyze_dir)

if __name__ == '__main__':
  if (len(sys.argv) < 3):
    print("Syntax: families_dir analyze_dir [method]")
    exit(1)
  print(" ".join(sys.argv))
  dataset_dir = sys.argv[1]
  analyze_dir = sys.argv[2]
  method = "JointSearch"
  if (len(sys.argv) > 3):
    method = sys.argv[3]
  analyse(dataset_dir, analyze_dir, method)




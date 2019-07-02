import os
import sys
import subprocess

sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/plotters')
import experiments as exp
import plot_two_arrays


def plot_evolution_ll(log_file):
  joint_ll = []
  reconciliation_ll = []
  sequences_ll = []
  iterations = []
  index = 0
  for line in open(log_file).readlines():
    if (not line.startswith("Likelihoods: ")):
      continue
    split = line.split(" ")
    joint_ll.append(float(split[3][:-1]))
    sequences_ll.append(float(split[6][:-1]))
    reconciliation_ll.append(float(split[9]))
    iterations.append(index)
    index += 1
  output = "likelihood_evolution.svg"
  plot_two_arrays.plot_two_arrays(output, iterations, reconciliation_ll, sequences_ll, "Search step", "Reconciliation LL", "Sequences LL", force_same_scale = True)

plot_evolution_ll("generax_logs.out")

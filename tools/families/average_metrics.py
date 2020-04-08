import sys
import os
sys.path.insert(0, os.path.join("tools", "print"))
sys.path.insert(0, os.path.join("script"))
sys.path.insert(0, os.path.join("tools", "print"))
import saved_metrics
import fam
import experiments as exp
from aligned_printer import AlignedPrinter


def print_average_metrics(metric_name, datadirs):
  metric_vector_dict = {}
  metric_vector_dict["generax-select-opt"] = []
  candidates = exp.generax_selector_candidates
  methods_in_best = open(candidates).read().replace("\n", "").lower().split(" ")
  print (methods_in_best)
  datasets_number = 0
  for datadir in datadirs:
    datadir_metrics = saved_metrics.get_metrics(datadir, metric_name)  
    if (datadir_metrics == None or len(datadir_metrics) == 0):
      continue
    datasets_number += 1
    bestMetric = 55555555
    for method in datadir_metrics:
      if (not method in metric_vector_dict):
        metric_vector_dict[method] = []
      val = float(datadir_metrics[method])
      metric_vector_dict[method].append(val)
      if (method.split(".")[0] in methods_in_best):
        bestMetric = min(bestMetric, val)
    for method in methods_in_best:
      if (method in datadir_metrics):
        bestMethod = min(bestMethod, float(datadir_metrics[method]))
    #metric_vector_dict["generax-select-opt"].append(bestMetric)
  if (datasets_number == 0):
    print("NO DATASET")
    sys.exit(1)
  print("Printing metric for " + str(datasets_number) + " datasets:")
  average_metrics = {}
  for method in metric_vector_dict:
    # compute the average only if all datasets contain the value
    v = metric_vector_dict[method]
    if (len(v) == datasets_number):
      average_metrics[method] = sum(v) / float(len(v))

  printer = AlignedPrinter()
  for method in average_metrics:
    printer.add(method + ":", str(average_metrics[method]))
  printer.sort_right_float()
  printer.display()

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print("Syntax: " + os.path.basename(__file__) + " metric_name directories")
    sys.exit(1)
  metric_name = sys.argv[1]
  datadirs = sys.argv[2:]
  print_average_metrics(metric_name, datadirs)



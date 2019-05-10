import os
import sys

import saved_metrics
sys.path.insert(0, os.path.join("tools", "plotters"))
import plot_two_arrays


ignore = ["generax-dtl-random", "generax-dtl-raxml", "generax-dl-random", "generax-weighted100"]

def compare_metrics(datadir, metric_name_1, metric_name_2, prefix):
  global ignore
  metric_dict_1 = saved_metrics.get_metrics(datadir, metric_name_1)
  metric_dict_2 = saved_metrics.get_metrics(datadir, metric_name_2)
  datax = []
  datax_labels = []
  alldata = []
  datay1 = []
  datay2 = []
  index = 0
  for method in metric_dict_1:
    if ((not method in metric_dict_2) or (method in ignore)):
      print("ignoring " + method)
      continue
    alldata.append((float(metric_dict_1[method]), float(metric_dict_2[method]), method))
  for data in reversed(sorted(alldata)):
    datay1.append(data[0])
    datay2.append(data[1])
    datax_labels.append(data[2])
    datax.append(index)
    index += 1
  output = prefix + "_" + metric_name_1.replace("_", "") + "_" + metric_name_2.replace("_", "") + ".png"
  plot_two_arrays.plot_two_arrays(output, datax, datay1, datay2, "method", metric_name_1, metric_name_2, xlabels = datax_labels)
  print("saved output in " + output)

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: datadir metric1 metric2 prefix")
    print("set prefix to show to show instead of saving")
    exit(1)

  datadir = sys.argv[1]
  metric_name_1 = sys.argv[2]
  metric_name_2 = sys.argv[3]
  prefix = sys.argv[4]
  compare_metrics(datadir, metric_name_1, metric_name_2, prefix)



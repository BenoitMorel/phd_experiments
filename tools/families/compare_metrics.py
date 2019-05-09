import os
import sys

import saved_metrics
sys.path.insert(0, os.path.join("tools", "plotters"))
import plot_two_arrays


def compare_metrics(datadir, metric_name_1, metric_name_2, output):
  metric_dict_1 = saved_metrics.get_metrics(datadir, metric_name_1)
  metric_dict_2 = saved_metrics.get_metrics(datadir, metric_name_2)
  datax = []
  datax_labels = []
  alldata = []
  datay1 = []
  datay2 = []
  index = 0
  for method in metric_dict_1:
    if (not method in metric_dict_2):
      continue
    alldata.append((float(metric_dict_1[method]), float(metric_dict_2[method]), method))
  for data in reversed(sorted(alldata)):
    datay1.append(data[0])
    datay2.append(data[1])
    datax_labels.append(data[2])
    datax.append(index)
    index += 1
  plot_two_arrays.plot_two_arrays(output, datax, datay1, datay2, "method", metric_name_1, metric_name_2, xlabels = datax_labels)

if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax: datadir metric1 metric2 output")
    print("set output to show to show instead of saving")
    exit(1)

  datadir = sys.argv[1]
  metric_name_1 = sys.argv[2]
  metric_name_2 = sys.argv[3]
  output = sys.argv[4]
  compare_metrics(datadir, metric_name_1, metric_name_2, output)



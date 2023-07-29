import os
import sys
import fast_rf_cells
sys.path.insert(0, os.path.join("tools", "plotters"))
from boxplot import BoxPlot


def get_run(paired_runs):
  return paired_runs.split("-")[1:]

def boxplot_rf(datadir, output, runs):
  rf_cells =  fast_rf_cells.load_rf_cells(datadir, False)
  boxplot = BoxPlot("Coucou", "Method", "RF distance")
  d = {} # d[run] = [values]
  for run in runs:
    d[run] = []
  for family in rf_cells:
    values = rf_cells[family]
    for run in runs:
      v = values["true.true - " + run]
      d[run].append(v[0] / v[1])
       
  for run in d:
    print(d[run])
    boxplot.add_elem(run.split("-")[0], d[run])
  boxplot.plot(output)
    


if __name__ == '__main__':
  if (len(sys.argv) < 3):
    print("Syntax: cells_datadir output run1 [run2 ...]")
    exit(1)
  datadir = sys.argv[1]
  output = sys.argv[2]
  runs = sys.argv[3:]
  boxplot_rf(datadir, output, runs)




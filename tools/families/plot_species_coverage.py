import os
import sys
sys.path.insert(0, 'tools/plotters')
sys.path.insert(0, 'tools/families')
import plot_histogram
import saved_metrics
import missing_data_estimation

def get_metrics(datadir):
  try:
    metrics = saved_metrics.get_metrics(datadir, "species_coverage")
    if (metrics == None):
      raise
    return metrics
  except:
    print("Failed to find species coverage information, computing it...")
    missing_data_estimation.estimate(datadir)
    return saved_metrics.get_metrics(datadir, "species_coverage")
  return None

def plot(datadir):
  output = os.path.basename(os.path.normpath(datadir)) + "_coverage.svg"
  metrics = get_metrics(datadir)
  x = []
  y = []
  for species in metrics:
    #x.append(species)
    x.append("")
    y.append(float(metrics[species]))
  plot_histogram.plot_histogram(x, y, sort_y = True, rotation = 90, output = output)
  print("Output file: " + output)

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  datadir = sys.argv[1]
  plot(datadir)






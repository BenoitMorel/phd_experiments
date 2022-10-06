import os
import sys
sys.path.insert(0, 'tools/plotters')
import plot_histogram
import glob



def get_likelihood(stats_file):
  line = open(stats).readlines()[0]



def plot_likelihood_distribution(results_dir, output):
  likelihoods = []
  xlabels = []
  for stats in glob.glob(os.path.join(results_dir, "*.stats")):
    #if ("14014_" in stats):
    #  continue
    line = open(stats).readlines()[0]
    sp = line.split()
    famll1 = float(sp[1])
    famll2 = float(sp[2])
    diff = famll2 - famll1
    if (diff > 0.0):
      likelihoods.append(diff)
      xlabels.append(len(xlabels))
  #data = {}
  #data["likelihoods"] = {}
  #for i in range(0, len(likelihoods)):
  #  data["likelihoods"][str(i)] = likelihoods[i]
  

  plot_histogram.plot_histogram(xlabels, likelihoods, output = output, xcaption = "MSAs", ycaption = "Likelihood diff", sort_y = True)







if __name__ == "__main__":
  if (len(sys.argv) < 2):
    print("syntax: python " + os.path.basename(__file__) + " results_dir output")
    sys.exit(1)
  results_dir = sys.argv[1]
  output = sys.argv[2]
  plot_likelihood_distribution(results_dir, output)
  




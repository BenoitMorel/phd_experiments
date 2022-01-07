import os
import sys
sys.path.insert(0, 'tools/plotters')
import plot_histogram
import random


def beta_reparametrized(mean, shape):
  alpha = mean * shape
  beta = (1 - mean) * shape
  return random.betavariate(alpha, beta)

def plot(mean, shape, points):
  output = "beta_" + str(mean) + "_" + str(shape) + "_" + str(points) + ".svg"
  x = []
  y = []
  for i in range(0, points):
    #x.append(species)
    x.append(i)
    y.append(beta_reparametrized(mean, shape))
  plot_histogram.plot_histogram(x, y, sort_y = True, rotation = 90, output = output)
  print("Output file: " + output)

if (__name__ == "__main__"):
  if (len(sys.argv) != 4):
    print("Syntax python " + os.path.basename(__file__) + " mean shape points")
    sys.exit(1)
  mean = float(sys.argv[1])
  shape = float(sys.argv[2])
  points = int(sys.argv[3])
  plot(mean, shape, points)







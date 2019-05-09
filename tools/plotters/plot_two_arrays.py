import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style("darkgrid")

def plot_two_arrays(output, datax, datay1, datay2, xlabel, ylabel1, ylabel2, xlabels = [], color1 = "red", color2 = "blue"):
  color1 = "tab:" + color1
  color2 = "tab:" + color2
  

  fig, ax1 = plt.subplots()
  ax1.set_xlabel("methods")
  ax1.set_ylabel(ylabel1, color=color1)
  ax1.plot(datax, datay1, color=color1)
  ax1.set_xticklabels(xlabels) 
  
  plt.xticks(rotation=45)
  plt.xticks(range(len(xlabels)), size='small')

  ax2 = ax1.twinx()
  ax2.set_ylabel(ylabel2, color=color2)
  ax2.plot(datax, datay2, color=color2)
  ax2.tick_params(axis="y", labelcolor = color2)
  

  fig.tight_layout()
  if (output == "show"):
    plt.show()
  else:
    plt.savefig(output)
    plt.close()

if (__name__ == "__main__"):
  t = ["meth1", "meth2", "meht3"]
  data1 = [1.0, 8.0, 1.0]
  data2 = [2.0, 1.0, 2.0]
  plot_two_arrays(t, data1, data2, "Method", "RF distance", "Adjacencies")




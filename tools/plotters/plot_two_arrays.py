import os
import sys
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set_style("darkgrid")


sns.set_style("darkgrid")

def plot_two_arrays(output, datax, datay1, datay2, xlabel, ylabel1, ylabel2, xlabels = [], color1 = "red", color2 = "blue", force_same_scale = False):
  color1 = "tab:" + color1
  color2 = "tab:" + color2
  

  fig, ax1 = plt.subplots()
  ax1.set_xlabel(xlabel)
  ax1.set_ylabel(ylabel1, color=color1)
  ax1.plot(datax, datay1, color=color1)
  ax1.set_xticklabels(xlabels) 
  ax1.tick_params(axis="y", labelcolor = color1)
    

  plt.xticks(rotation=45)
  plt.xticks(range(len(xlabels)), size='small')

  ax2 = ax1.twinx()
  ax2.set_ylabel(ylabel2, color=color2)
  ax2.plot(datax, datay2, color=color2)
  ax2.tick_params(axis="y", labelcolor = color2)
 
  if (force_same_scale):
    ylim1 = ax1.get_ylim()
    ylim2 = ax2.get_ylim()
    diff1 = ylim1[1] - ylim1[0]
    diff2 = ylim2[1] - ylim2[0]
    center1 = (ylim1[1] + ylim1[0]) / 2.0
    center2 = (ylim2[1] + ylim2[0]) / 2.0
    if (diff1 > diff2):
      ax2.set_ylim(center2 - diff1 / 2.0, center2 + diff1 / 2.0)
    else:
      ax1.set_ylim(center1 - diff2 / 2.0, center1 + diff2 / 2.0)
    
  fig.tight_layout()
  if (output == "show"):
    plt.show()
  else:
    print("saving figure in " + output)
    plt.savefig(output)
    plt.close()

if (__name__ == "__main__"):
  t = ["meth1", "meth2", "meht3"]
  data1 = [1.0, 8.0, 1.0]
  data2 = [2.0, 1.0, 2.0]
  plot_two_arrays("show", t, data1, data2, "Method", "RF distance", "Adjacencies")




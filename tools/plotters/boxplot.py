import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")
import pandas as pd






class BoxPlot:
  def __init__(self, title = None, ylabel = None):
    self._dict = {}
    self._title = title
    self._ylabel = ylabel

  def add_elem(self, xlabel, values):
    self._dict[xlabel] = values


  

  def plot(self, output):
    df = pd.DataFrame(data = self._dict)
    ax = sns.boxplot(data = df)
    plt.xticks(rotation = 45) 
    if (self._title != None):
      ax.set_title(self._title)
    if (self._ylabel != None):
      ax.set_ylabel(self._ylabel)
    ax.get_figure().tight_layout()
    if (output == "show"):
      plt.show()
    else:
      plt.savefig(output)
  

if (__name__ == "__main__"):
  
  np.random.seed(19680801)
  s = 10
  spread = np.random.rand(s) * 100
  center = np.ones(s) * 50
  flier_high = np.random.rand(s) * 100 + 100
  flier_low = np.random.rand(s) * -100
  bp = BoxPlot(title = "Example", ylabel = "ylabel")
  bp.add_elem("spread", spread)
  bp.add_elem("center", center)
  bp.add_elem("flier_high", flier_high)
  bp.add_elem("flier_low", flier_low)
  bp.plot("show")






import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")
import pandas as pd





class BoxPlot:
  def __init__(self, title = None, ylabel = None, order = None, max_y = -1.0):
    self._dict = {}
    self._title = title
    self._ylabel = ylabel
    self._order = order
    self._max_y = max_y

  def add_elem(self, xlabel, values):
    self._dict[xlabel] = values

  def plot(self, output):
    plt.clf()
    df = pd.DataFrame(data = self._dict)
    ax = sns.boxplot(data = df, order = self._order)
    if (self._max_y > 0.0):
      ax.set_ylim(0, self._max_y)
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

class GroupBoxPlot:
  """
  data is a dict of dict:   data[model][method] = values
  for each method, there is one boxplot per model with the values
  """
  def __init__(self, data, title = None, ylabel = None, xlabel = "x", hue_label = "hue", order = None):
    self._title = title
    self._xlabel = xlabel
    self._ylabel = ylabel
    self._hue_label = hue_label
    self._order = order
    self._x_number = 0
    hue_vector = []
    x_vector = []
    values_vector = []
    for hue in data:
      self._x_number = len(data[hue])
      for x in data[hue]:
        for value in data[hue][x]:
          hue_vector.append(hue)
          x_vector.append(x)
          values_vector.append(value)
    dataFrameDico = {}
    dataFrameDico[self._hue_label] = hue_vector
    dataFrameDico[self._xlabel] = x_vector
    dataFrameDico["values"] = values_vector
    self._df = pd.DataFrame(data = dataFrameDico)

  def plot(self, output):
    plt.clf()
    ax = sns.boxplot(data = self._df, hue = self._hue_label, x = self._xlabel, y = "values", order = self._order)
    if (self._xlabel == "x"):
      ax.set_xlabel('')    
    if (self._x_number > 4):
      print(self._x_number)
      plt.xticks(rotation = 45)
    if (self._title != None):
      ax.set_title(self._title)
    if (self._ylabel != None):
      ax.set_ylabel(self._ylabel)
    ax.get_figure().tight_layout()
    if (output == "show"):
      print("show")
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






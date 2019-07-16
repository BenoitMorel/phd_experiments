import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set_style("darkgrid")

def plot_histogram(xlabels, yvalues, title = None, xcaption = None, ycaption = None,start_at_min_y = False, reverse_bar = False, output = "show"):
    y_pos = np.arange(len(xlabels))
    fig, ax = plt.subplots()
    bottom = 0
    yvaluesbar = yvalues[:]
    if (reverse_bar):
      plop = 100000000
      bottom = -plop
      for i in range(0, len(yvalues)):
        yvaluesbar[i] = plop - abs(yvalues[i])
    plt.bar(y_pos, yvaluesbar, align='center', bottom = bottom)
    plt.xticks(y_pos, xlabels)
    plt.xticks(rotation=45)
    plt.xticks(range(len(xlabels)), size='small')
    if (xcaption != None):
      plt.xlabel(xcaption)
    if (ycaption != None): 
      plt.ylabel(ycaption)
    if (start_at_min_y):
      max_value = max(yvalues)
      min_value = min(yvalues)
      epsilon = (max_value - min_value) / 10.0
      ax.set_ylim(min_value - epsilon, max_value + epsilon)
    if (title != None):
      plt.title(title)
    fig.tight_layout()
    if (output == "show"):
        plt.show()
    else:
        print("Saving plot in " + output)
        plt.savefig(output)
"""
    data[category][class] = value
"""
def plot_grouped_histogram(data, title = None, xcaption = None, ycaption = None, cat_name = "categories", class_name = "classes", values_name = "values", kind = "bar", start_at_min_y = False, log_scale = False, x_order = None, output = "show"):
  values_vec = []
  categories_vec = []
  classes_vec = []
  min_value = 99999999999999
  max_value = -9999999999999
  for category in data:
    for classs in data[category]:
      value = data[category][classs]
      min_value = min(min_value, value)
      max_value = max(max_value, value)
      values_vec.append(value)
      categories_vec.append(category)
      classes_vec.append(classs)
  dataFrameDico = {}
  dataFrameDico[class_name] = classes_vec
  dataFrameDico[cat_name] = categories_vec
  dataFrameDico[values_name] = values_vec
  df = pd.DataFrame(data = dataFrameDico)
  f = sns.factorplot(data = df, x = class_name, hue = cat_name, y = values_name, kind = kind, order = x_order, legend = False)
  f.despine(left=True)
  #f._legend.set_title('')
  plt.gca().legend().set_title('')
  plt.gca().set_xlabel('')
  if (log_scale):
    f.set(yscale="log")
  plt.xticks(rotation = 90)
  if (start_at_min_y):
    epsilon = (max_value - min_value) / 10.0
    f.set(ylim=(min_value - epsilon, max_value + epsilon))
  if (title != None):
    f.fig.suptitle(title)
  f.fig.tight_layout()
  if (output == "show"):
    plt.show()
  else:
    print("Saving plot in " + output)
    plt.savefig(output)
  
    


if (__name__ == "__main__"):
    objects = ['Python', 'C++', 'Java', 'Perl', 'Scala', 'Lisp']
    performance = [10,8,6,4,2,1]
    #plot_histogram(objects, performance, "title",  "pif", "paf", "show")
    
    blabla = [5,4,3,2,2,1]
    data = {}
    data["blabla"] = {}
    data["performance"] = {}
    for i in range(0, len(objects)):
      data["blabla"][objects[i]] = blabla[i]
      data["performance"][objects[i]] = performance[i]
    plot_grouped_histogram(data = data, output = "show")



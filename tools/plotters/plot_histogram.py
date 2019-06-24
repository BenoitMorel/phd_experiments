import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set_style("darkgrid")

def plot_histogram(xlabels, yvalues, title, xcaption, ycaption, output):
    print(xlabels)
    print(yvalues)
    y_pos = np.arange(len(xlabels))
    fig, ax = plt.subplots()
    plt.bar(y_pos, yvalues, align='center')
    plt.xticks(y_pos, xlabels)
    plt.xticks(rotation=45)
    plt.xticks(range(len(xlabels)), size='small')
    plt.xlabel(xcaption)
    plt.ylabel(ycaption)
    plt.title(title)
    fig.tight_layout()
    if (output == "show"):
        plt.show()
    else:
        plt.savefig(output)
        plt.close()
"""
    data[category][class] = value
"""
def plot_grouped_histogram(title, xcaption = None, ycaption = None, output = "show"):
  values_vec = []
  categories_vec = []
  classes_vec = []
  for category in data:
    for classs in data[category]:
      values_vec.append(data[category][classs])
      categories_vec.append(category)
      classes_vec.append(classs)
  dataFrameDico = {}
  dataFrameDico["classes"] = classes_vec
  dataFrameDico["categories"] = categories_vec
  dataFrameDico["values"] = values_vec
  df = pd.DataFrame(data = dataFrameDico)
  #y_pos = np.arange(len(xlabels))
  sns.factorplot(data = df, x = "classes", hue = "categories", y = "values", kind = "bar")
  #plt.xticks(rotation = 45)
  #if (title != None):
  #  ax.set_title(title)
  #ax.get_figure().tight_layout()
  if (output == "show"):
    print("show")
    plt.show()
  else:
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
    plot_grouped_histogram(data, output = "show")



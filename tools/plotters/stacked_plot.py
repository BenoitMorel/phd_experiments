import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def stacked_plot(data, caption_title, stack_labels, xlabels, ylabel, log_scale = True, output = "show"):
  columns = [caption_title]
  columns.extend(xlabels)
  pd_data = []
  for i in range(0, len(data)):
    pd_row = [stack_labels[i]]
    for value in data[i]:
      pd_row.append(value)
    pd_data.append(pd_row)
  df = pd.DataFrame(columns = columns, data = pd_data)
  f = df.set_index(caption_title).T.plot(kind='bar', stacked=True)
  plt.ylabel(ylabel)
  if (log_scale):
    f.set(yscale="log")
  plt.xticks(rotation = 45)
    
  f.figure.tight_layout()
  if (output == "show"):
    plt.show()
  else:
    print("Saving plot in " + output)
    plt.savefig(output)
df = pd.DataFrame(columns=["yo","raxml-ng", "treercs", "generax"], 
                  data=[["tool",50, 10, 100],
                        ["precomputation",0,200, 0]])

#sns.set()
#df.set_index('yo').T.plot(kind='bar', stacked=True)
#plt.savefig("yop.svg")


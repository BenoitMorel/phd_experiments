import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

def plot_line(xvalues, yvalues_list, title, xcaption, ycaption, output, line_captions = None, sort_y = False, marker = "."):
    fig, ax = plt.subplots()
    lines_number = len(yvalues_list)
    assert(line_captions == None or len(line_captions) == lines_number)
    ax.tick_params(axis='both', which='major', labelsize=15)
    for i in range(0, lines_number):
        yvalues = yvalues_list[i]
        if (sort_y):
          hey = sorted(zip(yvalues, xvalues), reverse = True)
          tuples = zip(*hey)
          yvalues, values = [ list(tuple) for tuple in  tuples]
        if (line_captions != None):
          plt.plot(xvalues, yvalues, marker=marker, label = line_captions[i])
        else:
          plt.plot(xvalues, yvalues, marker=marker)
    fontsize = 20
    if (xcaption != None):
      plt.xlabel(xcaption, fontsize=fontsize)
    if (ycaption != None):
      plt.ylabel(ycaption, fontsize=fontsize)
    if (title != None):
      plt.title(title, fontsize=fontsize)
    plt.legend()
    fig.tight_layout()
    if (output == "show"):
        plt.show()
    else:
        plt.savefig(output)
        plt.close()

if (__name__ == "__main__"):
    xvalues = [8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0] 
    yvalues = [16975.0, 8690.0, 4752.0, 2658.0, 1561.0, 769.0, 495.0]
    cost_per_cores = yvalues[0] * xvalues[0]
    for i in range(0, len(xvalues)):
        yvalues[i] = cost_per_cores / yvalues[i]
    yvalues_list = [yvalues, xvalues]
    line_captions = ["GeneRax", "Theoretical optimum"]
    plot_line(xvalues, yvalues_list, "GeneRax parallel efficiency",  "Cores", "Speedup", "show", line_captions)


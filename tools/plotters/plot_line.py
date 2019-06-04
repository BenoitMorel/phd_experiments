import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

def plot_line(xvalues, yvalues_list, title, xcaption, ycaption, output, lines_captions = None):
    fig, ax = plt.subplots()
    lines_number = len(yvalues_list)
    assert(lines_captions == None or len(lines_captions) == lines_number)
    for i in range(0, lines_number):
        yvalues = yvalues_list[i]
        plt.plot(xvalues, yvalues, marker='.', label = lines_captions[i])
    plt.xlabel(xcaption)
    plt.ylabel(ycaption)
    plt.title(title)
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
    lines_captions = ["GeneRax", "Theoretical optimum"]
    plot_line(xvalues, yvalues_list, "GeneRax parallel efficiency",  "Cores", "Speedup", "show", lines_captions)


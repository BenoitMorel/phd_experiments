import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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


if (__name__ == "__main__"):
    objects = ['Python', 'C++', 'Java', 'Perl', 'Scala', 'Lisp']
    performance = [10,8,6,4,2,1]
    plot_histogram(objects, performance, "title",  "pif", "paf", "show")


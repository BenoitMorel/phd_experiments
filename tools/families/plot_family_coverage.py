import os
import sys
sys.path.insert(0, 'tools/plotters')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/mappings')
import plot_line
import plot_histogram
import saved_metrics
import ete3
import fam
import get_dico

def get_coverage(datadir):
  species_number = float(len(ete3.Tree(fam.get_species_tree(datadir)).get_leaves()))
  coverage = []
  for family in fam.get_families_list(datadir):
    d = get_dico.get_species_to_genes_family(datadir, family)
    coverage.append(float(len(d)) / species_number)
  return coverage

def plot(datadir):
  output = os.path.basename(os.path.normpath(datadir)) + "_family_coverage.svg"
  x = []
  y = get_coverage(datadir)
  i = 0
  for i in range(0, len(y)):
    x.append(i)
  title = None
  xcaption = None
  ycaption = None
  line_captions = None
  plot_line.plot_line(x, [y], title, xcaption, ycaption, output, line_captions, sort_y = True)# marker = None)
  #plot_histogram.plot_histogram(x, y, sort_y = True, rotation = 90, output = output)
  print("Output file: " + output)

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  datadir = sys.argv[1]
  plot(datadir)






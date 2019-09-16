import sys
sys.path.insert(0, 'tools/publications_generax/utils_plots')

import plot_scaling
import plot_rrf
import plot_boxplots
import plot_runtimes
import plot_ll



if (__name__ == "__main__"):
  plot_rrf.plot_simulated_metrics()
  plot_scaling.plot_scaling()  
  plot_boxplots.plot_model_boxplots()
  plot_runtimes.plot_runtimes()
  #plot_ll.plot_ll()


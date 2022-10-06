import os
import sys
sys.path.insert(0, 'tools/plotters')
import plot_line
import gather_likelihoods
  
def plot(rundir, totaltreenumber, treenumbers, output):
  
  yraxml = []
  ytreecombine = []
  for n in treenumbers:
    treecombine_run_dir = os.path.join(rundir, "treecombination" + str(totaltreenumber) + "_s" + str(n))
    results_dir = os.path.join(treecombine_run_dir, "results")
    llraxml, lltreecombine = gather_likelihoods.extract_likelihoods(results_dir, False)
    yraxml.append(llraxml)
    ytreecombine.append(lltreecombine)
  yvalues = [yraxml, ytreecombine]
  line_captions = ["RAxML-NG", "TreeCombine"]
  plot_line.plot_line(treenumbers, yvalues, "TreeCombine VS RAxML-NG",  "Starting trees", "Likelihood", output, line_captions)


if (__name__ == "__main__"):
  if (len(sys.argv) < 4):
    print("Syntax python " + os.path.basename(__file__) + " rundir totaltreenumber output [treenumber1, treenumber2, ...]")
    sys.exit(1)
  rundir = sys.argv[1]
  totaltreenumber = int(sys.argv[2])
  output = sys.argv[3]
  treenumbers = []
  for n in sys.argv[4:]:
    treenumbers.append(int(n))

  plot(rundir, totaltreenumber, treenumbers, output)

#if (__name__ == "__main__"):
#    xvalues = [10, 20, 40, 60, 100]
#    y1 = [-3467736.22, -3467673.3, -151166 + -3316446.07, -3219131.58 + -151166 + -97269.4, -3316103.73 + -151166]
#    y2 = [-3467600.08, -3467273.84, -151161 + -3316307.02, -3219043.11 + -151161 + -97236.7, -3315697.66 + -151161]
    
    #xvalues = [10, 20, 50, 75, 100]
    #y1 = [-3467740.39, -3467700.97, -3467588.53, -3467297.38, -3467269.73]
    #y2 = [-3467626.21, -3467558.15, -3467481.54, -3467184.22, -3467165.27]
#    yvalues = [y1, y2]
#    line_captions = ["RAxML-NG", "TreeCombine"]
#    plot_line.plot_line(xvalues, yvalues, "TreeCombine VS RAxML-NG",  "Starting trees", "Likelihood", "hey.svg", line_captions)




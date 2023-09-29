import os
import sys
sys.path.insert(0, os.path.join("tools", "ale"))
sys.path.insert(0, os.path.join("tools", "plots"))
import alerax_parser
import ale_parser
import boxplot


def plot_famrates(aledir, aleraxdir, output):
  data = {}
  aledata = {}
  print("Reading ALE rates...")
  alerates, aleevents = ale_parser.parse_all_rec_uml(aledir)
  print("Filling ALE data...")
  data["ALE"] = aledata
  aled = []
  alet = []
  alel = []
  aledata["D"] = aled
  aledata["L"] = alel
  aledata["T"] = alet
  for fam in alerates:
    v = alerates[fam]
    aled.append(v[0])
    alel.append(v[1])
    alet.append(v[2])
  print("Reading AleRax rates...")
  aleraxrates = alerax_parser.parse_per_family_rates(aleraxdir)
  print("Filling AleRax data...")
  aleraxdata = {}
  data["AleRax"] = aleraxdata
  aleraxd = []
  aleraxt = []
  aleraxl = []
  aleraxdata["D"] = aleraxd
  aleraxdata["L"] = aleraxl
  aleraxdata["T"] = aleraxt
  for fam in alerates:
    v = alerates[fam]
    aleraxd.append(v[0])
    aleraxl.append(v[1])
    aleraxt.append(v[2])
  print("Printing into " + output)
  plotter = boxplot.GroupBoxPlot(data, ylabel = "Rates", xlabel = "Events", violin = True, hue_label = "Method")
  plotter.plot(output)
  

if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " aledir aleraxdir output")
    sys.exit(0)
  aledir = sys.argv[1]
  aleraxdir = sys.argv[2]
  output = sys.argv[3]
  plot_famrates(aledir, aleraxdir, output)


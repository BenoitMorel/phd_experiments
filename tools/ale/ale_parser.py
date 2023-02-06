import os
import sys
sys.path.insert(0, os.path.join("tools", "plotters"))
from boxplot import BoxPlot
"""
 Parse files with a summary of the per-species events
 returns a dictionary: data[species] == [dup, trans, loss, orig, copy, num]
"""
def parse_events_summary(uml):
  lines = open(uml).readlines()
  events = {}
  for line in lines[1:]:
    sp = line.split()
    events[sp[0]] = [float(sp[1]), float(sp[2]), float(sp[3]), float(sp[4]), float(sp[5])]
  return events

def parse_rates_summary(uml):
  lines = open(uml).readlines()
  rates = {}
  for line in lines[1:]:
    sp = line.split()
    rates[sp[0]] = [float(sp[1]), float(sp[2]), float(sp[3])]
  return rates

def parse_rec_uml(uml):
  lines = open(uml).readlines()
  rates = []
  events = {}
  for line in lines:
    if (line.startswith("ML")):
      sp = line.split()
      for rate in sp[1:]:
        rates.append(float(rate))
      assert(len(rates) == 3)
      break
  for line in reversed(lines):
    if (line.startswith("#")):
      break
    sp = line.split()
    events[sp[1]] = [float(sp[2]), float(sp[3]), float(sp[4]), float(sp[5]), float(sp[6])] 
  return rates, events

def parse_all_rec_uml(uml_dir):
  all_events = {}
  all_rates = {}
  for fbase in os.listdir(uml_dir):
    if (fbase.startswith(".")):
      continue
    f = os.path.join(uml_dir, fbase)
    rates, events = parse_rec_uml(f)
    all_rates[fbase] = rates
    for species in events:
      if (not species in all_events):
        all_events[species]  = [0, 0, 0, 0, 0]
      for i in range(0, len(events[species])):
        all_events[species][i] += events[species][i]
  return all_rates, all_events

def check_almost_equal(list1, list2, epsilon = 0.0000001):
  assert(len(list1) == len(list2))
  for i in range(0, len(list1)):
    assert(abs(list1[i] - list2[i]) < epsilon)

def event_to_index(event):
  if (event == "Duplications"):
    return 0
  if (event == "Transfers"):
    return 1
  if (event == "Losses"):
    return 2

def check_tom(uml_dir, events_summary_file, rates_summary_file):
  events = parse_events_summary(events_summary_file)
  rates = parse_rates_summary(rates_summary_file)
  all_rates, all_events = parse_all_rec_uml(uml_dir)
  for species in events:
    check_almost_equal(events[species], all_events[species])
  for family in rates:
    check_almost_equal(rates[family], all_rates[family])


def plot_tom():
  datasets = []
  datasets.append(["Bacteria", "reports/supp_tom/bacteria/528_branchwise_events.txt", "reports/supp_tom/bacteria/528_rates.txt"])
  datasets.append(["Opisthokonta", "reports/supp_tom/fungi_metazoa/MLparams_branchwise_events.txt", "reports/supp_tom/fungi_metazoa/MLparams_rates.txt"])
  datasets.append(["Opisthokonta single copy", "reports/supp_tom/fungi_metazoa/sco_events.txt", "reports/supp_tom/fungi_metazoa/sco_rates.txt"])
  
  per_dataset_events = {}
  per_dataset_rates = {}
  per_dataset_verticality = {}
  per_dataset_td = {}
  for dataset_tuple in datasets:
    dataset = dataset_tuple[0]
    events = parse_events_summary(dataset_tuple[1])
    rates = parse_rates_summary(dataset_tuple[1])
    verticality = []
    td = []
    for species in events:
      d = events[species][0]
      t = events[species][1]
      s = events[species][4]
      v = (s + d) / (s + d + t)
      verticality.append(v)
      if (d > 0):
        td.append(t/d)
    per_dataset_events[dataset] = events
    per_dataset_rates[dataset] = rates
    per_dataset_verticality[dataset] = verticality
    per_dataset_td[dataset] = td


  
  print("Plotting events per branch:")
  for event in ["Duplications", "Losses", "Transfers"]:
    event_index = event_to_index(event)
    bp = BoxPlot(title = event, ylabel = "dataset", xlabel = "events per branch", max_x = 3000)
    for dataset_tuple in datasets:
      dataset = dataset_tuple[0]
      events = per_dataset_events[dataset]
      values = []
      for species in events:
        values.append(events[species][event_index])
      bp.add_elem(dataset, values)
    bp.plot("events_" + event + ".svg", True)
  

  print("Plotting rates per family:")
  for event in ["Duplications", "Losses", "Transfers"]:
    event_index = event_to_index(event)
    bp = BoxPlot(title = event, ylabel = "dataset", xlabel = "parameter per family", max_x = 1.50)
    for dataset_tuple in datasets:
      dataset = dataset_tuple[0]
      rates = per_dataset_rates[dataset]
      values = []
      for family in rates:
        values.append(rates[family][event_index])
      bp.add_elem(dataset, values)
    bp.plot("rates_" + event + ".svg", True)
  
  print("Plotting verticality:")
  bp = BoxPlot(title = event, ylabel = "dataset", xlabel = "verticality per branch", max_x = 1.0)
  for dataset_tuple in datasets:
    dataset = dataset_tuple[0]
    bp.add_elem(dataset, per_dataset_verticality[dataset])
  bp.plot("verticality.svg", True)
  
  print("Plotting td:")
  bp = BoxPlot(title = event, ylabel = "dataset", xlabel = "transfers/duplications per branch", max_x = 70.0)
  for dataset_tuple in datasets:
    dataset = dataset_tuple[0]
    bp.add_elem(dataset, per_dataset_td[dataset])
  bp.plot("td_ratio.svg", True)
  



if (__name__== "__main__"):
  max_args_number = 4
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + " uml_dir events_summary_file rates_summary_file")
    sys.exit(0)


  uml_dir = sys.argv[1]
  events_summary_file = sys.argv[2]
  rates_summary_file = sys.argv[3]

  #check_tom(uml_dir, events_summary_file, rates_summary_file)
  plot_tom()

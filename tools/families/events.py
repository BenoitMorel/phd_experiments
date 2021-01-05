import os
import sys
import pickle
import fam

"""
  all event counts aec: aec[model][generax_radius][event] = count
"""


def _save_all_event_counts(datadir, aec):
  out = fam.get_event_counts_file(datadir)
  with open(out, "w") as writer:
    for model in aec:
      for generax_radius in aec[model]:
        writer.write(model + "|")
        writer.write(str(generax_radius) + "|")
        event_counts = aec[model][generax_radius]
        l = []
        for event in event_counts:
          l.append(event + "=" + str(event_counts[event]))
        writer.write(",".join(l))
        writer.write("\n")

def _add_event_counts(aec, model, generax_radius, event_counts):
  if (not model in aec):
    aec[model] = {}
  aec[model][generax_radius] = event_counts

def get_all_event_counts(datadir):
  f = fam.get_event_counts_file(datadir) 
  aec = {}
  if (not os.path.isfile(f)):
    print("Warning: " + f + " is not a file")
    return aec
  for line in open(f).readlines():
    sp = line.split("|")
    if (len(sp) < 3):
      continue
    model = sp[0]
    generax_radius = int(sp[1])
    event_counts = {}
    for p in sp[2].split(","):
      sp2 = p.split("=")
      event_counts[sp2[0]] = int(sp2[1])
    _add_event_counts(aec, model, generax_radius, event_counts)
  return aec

"""
  event_counts is a dictionnary event_label->count
  for instance event_counts["T"] == 5 for 5 transfers
"""
def update_event_counts(datadir, model, generax_radius, event_counts):
  aec = get_all_event_counts(datadir) 
  _add_event_counts(aec, model, generax_radius, event_counts)
  _save_all_event_counts(datadir, aec)

def print_event_counts(datadir):
  aec = get_all_event_counts(datadir)
  if (aec == {}):
    print("No event counts")
    return
  for model in aec:
    for generax_radius in aec[model]:
      toprint = "counts for model " + model 
      toprint += " and radius " + str(generax_radius)
      toprint += "\n"
      event_counts = aec[model][generax_radius]
      for event in event_counts:
        toprint += "  " + event + "\t= " + str(event_counts[event]) + "\n"
      print(toprint)

if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " datadir")
    sys.exit(1)
  datadir = sys.argv[1]
  print_event_counts(datadir)



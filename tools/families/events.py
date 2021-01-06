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

def get_D_S(event_counts):
  S = event_counts["S"] + event_counts["SL"]
  D = event_counts["D"]
  return float(D) / float(S)

def get_D_lf(event_counts):
  D = event_counts["D"]
  leaves = event_counts["Leaf"]
  return float(D) / float(leaves)

def get_T_S(event_counts):
  S = event_counts["S"] + event_counts["SL"]
  T = event_counts["T"] + event_counts["TL"]
  return float(T) / float(S)

def get_T_lf(event_counts):
  T = event_counts["T"] + event_counts["TL"]
  leaves = event_counts["Leaf"]
  return float(T) / float(leaves)

def get_D_L(event_counts):
  L = event_counts["TL"] + event_counts["SL"]
  D = event_counts["D"]
  return float(D)/float(L)

def print_event_counts(aec):
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
      toprint += "  D/S\t= " + str(get_D_S(event_counts)) + "\n"
      toprint += "  D/lf\t= " + str(get_D_lf(event_counts)) + "\n"
      toprint += "  D/L\t= " + str(get_D_L(event_counts)) + "\n"
      toprint += "  T/S\t= " + str(get_T_S(event_counts)) + "\n"
      toprint += "  T/lf\t= " + str(get_T_lf(event_counts)) + "\n"
      print(toprint)

def print_event_counts_datadir(datadir):
  print_event_counts(get_all_event_counts(datadir))

def sum_event_counts(datadirs):
  sum_aec = get_all_event_counts(datadirs[0])
  for datadir in datadirs[1:]:
    aec = get_all_event_counts(datadir)
    if (aec == {}):
      continue
    assert(aec.keys() == sum_aec.keys())
    for model in aec:
      assert(aec[model].keys() == sum_aec[model].keys())
      for generax_radius in aec[model]:
        ec = aec[model][generax_radius]
        sum_ec = sum_aec[model][generax_radius]        
        for event in sum_ec:
          sum_ec[event] += ec[event]
  print_event_counts(sum_aec)

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + " datadirs")
    sys.exit(1)
  if (len(sys.argv) == 2):
    datadir = sys.argv[1]
    print_event_counts_datadir(datadir)
  else:
    datadirs = sys.argv[1:]
    sum_event_counts(datadirs)



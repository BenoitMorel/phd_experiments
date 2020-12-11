import sys
import os
import glob



def extract(generaxdir):
  events = {}
  events["S"] = 0
  events["SL"] = 0
  events["D"] = 0
  events["T"] = 0
  events["TL"] = 0
  events["L"] = 0
  events["Leaf"] = 0
  pattern = os.path.join(generaxdir, "reconciliations", "*_eventCounts.txt")
  event_count_files = glob.glob(pattern)
  for f in event_count_files:
    lines = open(f).readlines()
    for line in lines:
      sp = line.replace("\n", "").split(":")
      events[sp[0]] = events[sp[0]] + int(sp[1])
  return events

def run(generaxdir):
  events = extract(generaxdir)
  for event in events:
    print(event + ":" + str(events[event]))


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax python " + os.path.basename(__file__) + " generaxdir")
    sys.exit(1)
  run(sys.argv[1])

import os
import sys
sys.path.insert(0, 'scripts')
import experiments as exp
import random
import shutil



def sample_from_grove_aux(model, size, sampled):
  d = os.path.join(exp.raw_datasets_root, "grove_" + str(size) + "_" + model)
  datasets = os.listdir(d)
  log = os.path.join(d, sampled, "log_0.txt")
  lines = open(log).readlines()
  states = set(["A", "C", "G", "T"])
  symrates = {}
  q = {}
  for state in states: 
    symrates[state] = {}
    symrates[state][state] = 1.0
  freqs = {}
  for line in lines:
    if (line.startswith("freq")):
      sp = line.split()
      st = sp[1][-3]
      freq = float(sp[-1])
      freqs[st] = freq
    if (line.startswith("rate")):
      sp = line.split()
      st1 = sp[1]
      st2 = sp[3][0]
      rate = float(sp[-1])
      symrates[st1][st2] = rate
      symrates[st2][st1] = rate
  if (len(freqs) == 0):
    shutil.rmtree(os.path.join(d,sampled))
    print("invalid model " + sampled)
    return sample_from_grove(model, size)
  freqsum = 0.0
  for state in states:
    freqsum += float(freqs[state])
  for state in states:
    freqs[state] = float(freqs[state]) / freqsum
  for st1 in states:  
    for st2 in states:
      q[st1 + st2] = symrates[st1][st2] * freqs[st2]
  print(freqs)
  return (q, symrates, freqs)

def sample_from_grove(model, size):
  d = os.path.join(exp.raw_datasets_root, "grove_" + str(size) + "_" + model)
  datasets = os.listdir(d)
  sampled = random.choice(datasets)
  return sample_from_grove_aux(model, size, sampled)

def clean_grove(model, size):
  d = os.path.join(exp.raw_datasets_root, "grove_" + str(size) + "_" + model)
  datasets = os.listdir(d)
  for dataset in datasets:
    try:
      res = sample_from_grove_aux(model, size, dataset)
    except:
      print("error in " + os.path.join(d, dataset))
      shutil.rmtree(os.path.join(d, dataset))

def sample_from_grove_safe(model, size):
  try:
    return sample_from_grove(model, size)
  except:
    print("Error, resampling")
    return sample_from_grove_safe(model, size)

if (__name__ == "__main__"):
  if (len(sys.argv) < 1):
    print("Syntax python " + os.path.basename(__file__))
    sys.exit(1)
  #clean_grove("gtr", 100)
  sample_from_grove_safe("gtr", 100)


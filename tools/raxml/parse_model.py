import os
import sys


def extract_model(raxml_model_file):
  res = lambda:0
  line = open(raxml_model_file).readlines()[0]
  split = line.split("+")
  str1 = split[0]
  str2 = split[1]
  res.model = str1[0:str1.find("{")]
  res.rates = str1[str1.find("{")+1:str1.find("}")].split("/")
  res.frequencies = str2[str2.find("{")+1:str2.find("}")].split("/")
  return res

def model_to_seqgen_cmd(model):
  cmd = []
  cmd.append("-m")
  cmd.append(model.model)
  cmd.append("-r")
  cmd.extend([str(i) for i in model.rates])
  cmd.append("-f")
  cmd.extend([str(i) for i in model.frequencies])
  return cmd


if (__name__ == "__main__"):
  if (len(sys.argv) != 2):
    print("Syntax: python parse_model.py model_file")
    exit(1)
  model = extract_model(sys.argv[1])
  print("model: " + str(model.model))
  print("rates: " + str(model.rates))
  print("frequencies: " + str(model.frequencies))

  seqgen_cmd = model_to_seqgen_cmd(model)
  print(" ".join(seqgen_cmd))


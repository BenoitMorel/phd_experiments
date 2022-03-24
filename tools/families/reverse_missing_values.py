import os
import sys
import shutil
import re
sys.path.insert(0, 'scripts')
import experiments as exp
import fam_data
import fam



def reverse(datadir):
  print(datadir)
  params = set(["ms", "mf"])
  sp = datadir.split("_")
  newsp = []
  for i in range(0, len(sp)):
    elem = sp[i]
    param_name =  re.sub("[0-9]*[\.]*[0-9]*", "", elem)
    param_name = re.sub("e-", "", param_name)
    if (param_name in params):
      value =  float(elem[len(param_name):])
      new_value = 1.0 - value
      elem = param_name + str(new_value)
    newsp.append(elem)
  new_datadir = "_".join(newsp)
  print(new_datadir)
  shutil.move(datadir, new_datadir)
    
def reverse_datadirs(datadirs):
  for datadir in datadirs:
    reverse(datadir)

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + " datadirs")
    sys.exit(1)
  datadirs = sys.argv[1:]
  reverse_datadirs(datadirs)



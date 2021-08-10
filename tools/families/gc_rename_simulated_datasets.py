import os
import sys
import shutil

def add_gc(datadir):
  print(datadir)
  if ("gc" in datadir):
    return
  key = "_p0"
  sp = datadir.split(key)
  newdatadir = sp[0] + "_gc0.0" + key + sp[1]
  print(newdatadir)
  shutil.move(datadir, newdatadir)
  
def rm_gc(datadir):
  print(datadir)
  key = "_gc0.0"
  if (not key in datadir):
    return
  sp = datadir.split(key)
  assert(len(sp) > 1)
  newdatadir = sp[0] + sp[1]
  print(newdatadir)
  shutil.move(datadir, newdatadir)

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + " datasets")
    sys.exit(1)
  for datadir in sys.argv[1:]:
    add_gc(datadir)




import sys
import os

def mysplit(lines):
  split = []
  block = []
  for line in lines:
    if (len(line.strip()) == 0):
      split.append(block)
      block = []
    else:
      block.append(line)
  return split

def is_ok(block, ignore):
  for line in block:
    for value in ignore:
      if value in line:
        return False
  return True

def remove_blocks(blocks, ignore):
  new_blocks = []
  for block in blocks:
    if (is_ok(block, ignore)):
      new_blocks.append(block)
  return new_blocks

def analyze(valgrind_file, ignore):
  lines = open(valgrind_file).readlines()
  lines = [s.split("==")[2][:-1] for s in lines if len(s.split("==")) == 3]
    
  blocks = mysplit(lines)
  blocks = remove_blocks(blocks, ignore)
  for b in blocks:
    for line in b:
      print(line)
    print("-----------------------")
  #print(lines)

if (__name__ == "__main__"): 
  if (len(sys.argv) != 2): 
    print("Syntax: python " + os.path.basename(__file__) + " valgrind_file")
    exit(1)

  ignore = ["PMPI_Init", "start_thread", "ompi_mpi_init", "PMPI_Comm_split", "PMPI_Comm_dup", "dlopen", "mpiSplit", "ParallelContext::init", "opal_db_base_fetch", "doWork"]

  valgrind_file = sys.argv[1]
  analyze(valgrind_file, ignore)





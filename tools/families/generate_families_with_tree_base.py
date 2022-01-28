import os
import sys
sys.path.insert(0, 'tools/families')
import generate_families_with_partitionned_msa as generate
import fam




if (__name__ == "__main__"): 
  if (len(sys.argv) < 2): 
    print("Syntax: python " + os.path.basename(__file__) + " treebaseIDs")
    exit(1)
  IDs = sys.argv[1:]
  for ID in IDs:
    treebase_dir = '/hits/fast/cme/TreeBase/'
    msa_file = os.path.join(treebase_dir, ID + ".phy")
    partition_file = os.path.join(treebase_dir, ID + ".part")
    is_dna = True
    species_tree_file = "NONE"
    datadir = fam.get_datadir("treebase_" + ID)
    failures = 0
    try:
      generate.generate(msa_file, partition_file, is_dna, species_tree_file, datadir)
    except:
      print("Failed to analyze " + ID)
      failures = failures + 1
    print("NUMBER OF FAILURES: " + str(failures))


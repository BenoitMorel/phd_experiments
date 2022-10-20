import os
import sys
sys.path.insert(0, 'tools/trees')
import read_tree
import relabel_from_tree 

 
def count_expected_average_compatible(d):
  ok = 0
  ko = 0
  for entry1 in d:
    for entry2 in d:
      from_start = d[entry1][0]
      to_end = d[entry2][1]
      compatible = from_start < to_end
      if (compatible):
        ok += 1
      else:
        ko += 1
  return float(ok) / float(ok + ko)


def count_transfers(generax_dir, species_tree_path):
  generax_tree = os.path.join(generax_dir, "species_trees", "inferred_species_tree.newick")
  tree = relabel_from_tree.load_relabeled(species_tree_path, generax_tree)
  root = tree
  index = 0
  d = {}
  for node in tree.traverse("postorder"):
    depth = node.get_distance(root)
    d[node.name] = (depth - node.dist, depth)

  print("expected average compatible: " + str(count_expected_average_compatible(d)))
  transfers_path = os.path.join(generax_dir, "per_species_pair_transfers.txt")
  total_transfers = 0
  total_ok_transfers = 0
  total_ko_transfers = 0
  total_contemporary_transfers = 0
  for line in open(transfers_path).readlines():
    sp = line.split()
    from_sp = sp[0]
    to_sp = sp[1]
    count = int(sp[2])
    total_transfers += count
    from_start = d[from_sp][0]
    from_end = d[from_sp][1]
    to_start = d[to_sp][0]
    to_end = d[to_sp][1]
    compatible = from_start < to_end
    
    plot_threshold = 500
    if (compatible):
      total_ok_transfers += count
      if (count > plot_threshold):
        print("OK " + line[:-1])
    else:
      if (count > plot_threshold):
        print("KO " + line[:-1])
      total_ko_transfers += count
    
    contemporary = compatible and to_start < from_end
    if (contemporary):
      total_contemporary_transfers += count


  proportion_ok = float(total_ok_transfers) / float(total_transfers)
  print("Total number of transfers: " + str(total_transfers))
  print("Total number of contemporary transfers: " + str(total_contemporary_transfers))
  print("Total number of OK transfers: " + str(total_ok_transfers) + " (" + str(proportion_ok) + ")")
  print("Total number of KO transfers: " + str(total_ko_transfers))
  

if (__name__ == "__main__"):
  if (len(sys.argv) < 2):
    print("Syntax python " + os.path.basename(__file__) + "generax_dir species_tree_path")
    sys.exit(1)
  generax_dir = sys.argv[1]
  species_tree_path = sys.argv[2]
  count_transfers(generax_dir, species_tree_path)


import os
import sys
import ete3




def rename(input_tree, output_tree):
  d = {}
  d["Bambusicola_thoracicus"] = "bambusicola_thoracicus"
  d["Caloperdix_oculeus"] = "caloperdix_oculea"
  d["Excalfactoria_chinensis"] = "excalfactoria_chinensis"
  d["Falcipennis_canadensis"] = "falcipennis_canadensis"
  d["arborophila_torqueola"] = "arborophila_torqeola"
  d["catreus_wallichii"] = "catreus_wallichi"
  d["melanoperdix_niger"] = "melanoperdix_nigra"
  d["meleagris_ocellata"] = "meleagris_ocellatus"
  d["perdix_dauurica"] = "perdix_dauuricae"
  d["pternistis_ahantensis"] = "pternistis_ahatensis"
  tree = ete3.Tree(input_tree, format = 1)
  for leaf in tree.get_leaves():
    if (leaf.name in d):
      leaf.name = d[leaf.name]
  tree.write(format = 1, outfile = output_tree)

if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax python gallimissing_rename_tree.py input_tree output_tree")
    sys.exit(1)
  input_tree = sys.argv[1]
  output_tree = sys.argv[2]
  rename(input_tree, output_tree)

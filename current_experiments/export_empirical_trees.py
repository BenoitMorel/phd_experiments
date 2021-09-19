import os
import sys
sys.path.insert(0, os.path.join("tools", "families"))
import fam
import translate_species_tree

methods = []
methods.append(("astralpro", "astralpro"))
methods.append(("generax-MiniNJ-prune-fam", "speciesraxprune"))
methods.append(("generax-MiniNJ-fam", "speciesrax"))
methods.append(("njrax-MiniNJ", "mininj"))
methods.append(("duptree", "duptree"))
methods.append(("fastmulrfs-single", "fastmulrfs"))


triplets = []
triplets.append(("pdb_plants23", "plants23", "bestAA"))
triplets.append(("pdb_fungi59", "fungi59", "bestAA"))
triplets.append(("pdb_fungi60", "fungi60", "bestAA"))
triplets.append(("pdb_vertebrates21", "vertebrates21", "true"))
triplets.append(("apro_plants", "plants83", "true"))
triplets.append(("apro_fungi", "fungi16", "true"))
triplets.append(("vertebrates188", "vertebrates188", "GTR+G"))
triplets.append(("aa_ensembl_98_ncrna_primates", "primates13", "GTR+G"))
triplets.append(("cyano_empirical", "cyanobacteria36", "LG+G+I"))


if (__name__ == "__main__"):
  
  outputdir = "empirical_species_trees" 
  
  os.mkdir(outputdir)
  for triplet in triplets:
    dataset, name, model = triplet
    datadir = fam.get_datadir(dataset)
    tripletdir = os.path.join(outputdir, name)
    os.mkdir(tripletdir)
    for method in methods:
      suffix = "_raxml-ng"
      if (model == "true"):
        suffix = "_true"
      species_tree = fam.get_species_tree(datadir, model, method[0] + suffix)
      out_species_tree = os.path.join(tripletdir, method[1] + ".newick")
      print(species_tree)
      translate_species_tree.dump_into(species_tree, out_species_tree) 

  

import os
import sys
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import fam_data
import fam
import experiments as exp
import saved_metrics
import rf_distance
from ete3 import Tree

# list of tuples (long method name, short method name)
methods_tuples = []
methods_tuples.append(("generax-MiniNJ-prune-fam", "S"))
methods_tuples.append(("njrax-MiniNJ", "M"))
methods_tuples.append(("astralpro", "A"))
methods_tuples.append(("fastmulrfs-single", "F"))
#methods_tuples.append(("duptree", "D"))

methods = []
for m in methods_tuples: 
  methods.append(m[1])

# list of tuples (dataset, paper_name, subst model, gene tree method)
datasets_tuples = []
#datasets_tuples.append(("aa_ensembl_98_ncrna_primates", "primates13", "GTR+G", "raxml-ng"))
#datasets_tuples.append(("cyano_empirical", "cyanobacteria36", "LG+G+I", "raxml-ng"))
#datasets_tuples.append(("pdb_vertebrates21", "vertebrates21", "true", "true"))
#datasets_tuples.append(("pdb_fungi59", "fungi59", "bestAA", "raxml-ng"))
#datasets_tuples.append(("pdb_fungi60", "fungi60", "bestAA", "raxml-ng"))
#datasets_tuples.append(("aa_plants", "plants23", "bestAA", "raxml-ng"))
#datasets_tuples.append(("apro_plants", "plants83", "true", "true"))
#datasets_tuples.append(("vertebrates188", "vertebrates188", "GTR+G", "raxml-ng"))
datasets_tuples.append(("tom", "tom", "true", "true"))

tex = "\\begin{table*}\\begin{center}\n"
tables_per_line = 2
index = 0
for dataset in datasets_tuples:
  
  trees = {}
  datadir = fam.get_datadir(dataset[0])
  caption = dataset[1]
  tex += "\\subfloat[" + caption + "]{\n"
  tex += "\\begin{tabular}{|"
  for m in methods_tuples:

    
    species_method = m[0] + "_" + dataset[3]
    tree_path = fam.get_species_tree(datadir, dataset[2], species_method)
    try:
      trees[m[1]] = Tree(tree_path, format = 1)
    except:
      print("Cannot read tree " + tree_path)
      sys.exit(1)
  for m in methods_tuples:
    tex += "c|"
  tex += "c|"
  tex += "}\n\\hline " + " & " + " & ".join(methods) + " \\\\\n"
  for m1 in methods_tuples:
    tree1 = trees[m1[1]]
    distances = []
    for m2 in methods_tuples:
      tree2 = trees[m2[1]]
      rrf = rf_distance.get_rf(tree1, tree2)
      distances.append(str(rrf)[:5])
    tex += "\\hline " + m1[1] + " & " + " & ".join(distances) +  " \\\\\n"
  tex += "\\hline\n"
  tex += "\\end{tabular}\n"
  tex += "\label{table:" + dataset[1] + "} \n}\n"
  if (index % tables_per_line == 1):
    tex += "\n"
  index += 1

tex += "\end{center} \end{table*}\n" 
print(tex)

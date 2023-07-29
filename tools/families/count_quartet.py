import os
import sys
import subprocess
import shutil
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
sys.path.insert(0, 'tools/mappings')
import experiments as exp
import fam
import get_dico


def count(datadir, method, model, A, B, C, D):
  frequencies = [0, 0, 0 , 0]
  for family in fam.get_families_list(datadir):
    species_to_genes = get_dico.get_species_to_genes_family(datadir, family)
    
    try:
      genesA = species_to_genes[A]
      genesB = species_to_genes[B]
      genesC = species_to_genes[C]
      genesD = species_to_genes[D]
      ok = True
      if (len(genesA) != 1 or len(genesB) != 1 or len(genesC) != 1 or len(genesD) != 1):
        continue
      tree_path = fam.get_gene_tree_path(datadir, family, method, model)
      command = []
      command.append(exp.quartet_counter_exec)
      command.append(tree_path)
      command.append(genesA[0])
      command.append(genesB[0])
      command.append(genesC[0])
      command.append(genesD[0])
      res = subprocess.check_output(command)
      frequencies[int(res)] += 1
      print(frequencies)
    except:
      print("skip " + family)
      continue
  s = frequencies[1] + frequencies[2] + frequencies[3]
  print(str(float(frequencies[1] / float(s))))
  print("((" + A + ", " + B + "), (" + C + "," + D + "));")
  print(str(float(frequencies[2] / float(s))))
  print("((" + A + ", " + C + "), (" + B + "," + D + "));")
  print(str(float(frequencies[3] / float(s))))
  print("((" + A + ", " + D + "), (" + B + "," + C + "));")

if (__name__ == "__main__"):
  if (len(sys.argv) < 8):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_method subst_model A B C D")
    sys.exit(1)
  datadir = sys.argv[1]
  method = sys.argv[2]
  model = sys.argv[3]
  A = sys.argv[4]
  B = sys.argv[5]
  C = sys.argv[6]
  D = sys.argv[7]
  count(datadir, method, model, A, B, C, D)



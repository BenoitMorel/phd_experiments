import os
import sys
import shutil
sys.path.insert(0, os.path.join("tools", "families"))
import fam

def export_note(output_dir):
  with open(os.path.join(output_dir, "README.txt"), "w") as writer:
    writer.write("This directory contains the datasets used in GeneRax paper.\n")
    writer.write("The cyanobacteria datasets comes from the ALE paper.\n")
    writer.write("We extracted the primates dataset from ENSEMBL.\n")
    writer.write("We generated the jsim datasets with jprime and seqgen.\n")
    writer.write("\n")
    writer.write("Each directory follows the exact same structure.\n")
    writer.write("- species_tree contains the species trees.\n")
    writer.write("- families contains the gene families.\n")
    writer.write("- alignments contains symlink to the alignments (the real files are in the family directories).\n")
    writer.write("\n")
    writer.write("Each gene family directory contains:\n")
    writer.write("- the MSA alignment (either fasta or phylip).\n")
    writer.write("- the gene to species mapping file.\n")
    writer.write("- gene_trees: a directory with all the trees we inferred for the paper, including the true tree for simulated datasets (for empirial datasets, the \"true\" trees are the trees from the database and should not be seen as the ground truth!!).\n")
    writer.write(".\n")

def export_dataset(dataset, output_dir):
  datadir = fam.get_datadir(dataset)
  newdatadir = os.path.join(output_dir, dataset)
  os.mkdir(newdatadir)
  ignore = shutil.ignore_patterns("misc")
  print("  copy species tree")
  shutil.copytree(fam.get_species_dir(datadir), fam.get_species_dir(newdatadir))
  print("  copy alignments directory")
  shutil.copytree(fam.get_alignments_dir(datadir), fam.get_alignments_dir(newdatadir), symlinks = True)
  print("  copy families directory")
  shutil.copytree(fam.get_families_dir(datadir), fam.get_families_dir(newdatadir), ignore=ignore)

def export_all_data(output_dir):
  print("Removing previous files...")
  try:
    shutil.rmtree(output_dir)
  except:
    pass
  os.mkdir(output_dir)
  export_note(output_dir)
  datasets = []
  for dataset in os.listdir(os.path.join(fam.get_datasets_family_path())):
    if (dataset.startswith("jsim")):
      datasets.append(dataset)
  datasets.append("cyano_empirical")
  datasets.append("cyano_simulated")
  datasets.append("ensembl_96_ncrna_primates")
  for dataset in datasets:
    print("Exporting dataset " + dataset)
    export_dataset(dataset, output_dir)

if (__name__ == "__main__"):
  output_dir = "extracted_data"
  export_all_data(output_dir)
  


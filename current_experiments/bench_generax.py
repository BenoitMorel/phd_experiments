import subprocess


def run_generax(dataset, species_tree, do_prune, per_fam, constrained, seed, rec_model = "UndatedDTL"):
  command = []
  command.append("python")
  command.append("scripts/generax/launch_speciesrax.py")
  command.append(dataset)
  command.append("GTR")
  command.append(species_tree)
  command.append("fasttree") # gene trees
  command.append("normald")
  command.append("45") # cores
  if (do_prune):
    command.append("--prune-species-tree")
  if (per_fam):
    command.append("--per-family-rates")
  if (not constrained):
    command.append("--unconstrained-species-search")
  command.append("--rec-model")
  command.append("UndatedDTL")
  command.append("--do-not-reconcile")
  command.append("--seed")
  command.append(str(seed))
  subprocess.check_call(command)

datasets = []
datasets.append("aa_ensembl_98_ncrna_primates")
#datasets.append("aa_ensembl_98_ncrna_sauropsids")
datasets.append("aa_ensembl_98_ncrna_lowprimates")
datasets.append("aa_ensembl_98_ncrna_mammals")
datasets.append("aa_ensembl_98_ncrna_vertebrates")
#datasets.append("aa_ensembl_98_ncrna_fishes")


for dataset in datasets:
  if (True):
    run_generax(dataset, "MiniNJ", True, True, False, 1)
    run_generax(dataset, "random", True, True, False, 1)
    run_generax(dataset, "random", True, True, False, 2)
    run_generax(dataset, "random", True, True, False, 3)
    




library(ape)
library(phangorn)

dist_family <- function(family_name, ale_rf)
{
  prefix <- paste("/hits/basement/cme/morel/github/BenoitDatasets/families/cyano_simulated/families", family_name, sep="/")
  prefix <- paste(prefix, "gene_trees", sep="/")
  true_tree <- read.tree(paste(prefix, "true.true.geneTree.newick", sep="/"))
  generax_tree <- read.tree(paste(prefix, "generax-dtl-random.LG+G+I.geneTree.newick", sep="/"))
  ale_tree <- read.tree(paste(prefix, "ale-dtl.LG+G+I.geneTree.newick", sep="/"))
  print(RF.dist(true_tree, generax_tree))
  print(KF.dist(true_tree, generax_tree))
  print(mean(RF.dist(true_tree, ale_tree)))
  print(mean(KF.dist(true_tree, ale_tree)))
  print("")
  run_rf[0] <- run_rf[0] + mean(RF.dist(true_tree, ale_tree))
}


families_path <- "/hits/basement/cme/morel/github/BenoitDatasets/families/cyano_simulated/families"
families <- list.files(families_path)

true_tree_str <- "true.true.geneTree.newick"
runs <- c("generax-dtl-random.LG+G+I.geneTree.newick", 
          "treerecs.LG+G+I.geneTree.newick",
          "phyldog.LG+G+I.geneTree.newick",
          "eccetera.LG+G+I.geneTree.newick",
          "notung90.LG+G+I.geneTree.newick",
          "raxml-ng.LG+G+I.geneTree.newick",
          "ale-dtl.LG+G+I.geneTree.newick")


for (run in runs) {
  run_rf = 0.0
  run_kf = 0.0
  for (family in families) {
    prefix <- paste(families_path, family, sep="/")
    prefix <- paste(prefix, "gene_trees", sep="/")
    true_tree <- unroot(read.tree(paste(prefix, true_tree_str, sep="/")))
    generax_tree <- unroot(read.tree(paste(prefix, run, sep="/")))
    run_rf <- run_rf + mean(RF.dist(true_tree, generax_tree, normalize = TRUE))
    run_kf <- run_kf + mean(KF.dist(true_tree, generax_tree))
  }
  average_run_rf = run_rf / length(families)
  average_run_kf = run_kf / length(families)
  print(run)
  print(paste("average RF:", average_run_rf))
  print(paste("average KF:", average_run_kf))
  print("")
}




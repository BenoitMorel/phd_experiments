
main <- function(tree1_path, tree2_path)
{
  tree1 <- read.tree(tree1_path)
  tree2 <- read.tree(tree2_path)
  run_kf <- mean(KF.dist(tree1, tree2))
  print(run_kf)
}

args = commandArgs(trailingOnly=TRUE)
tree1_path <- args[1]
tree2_path <- args[2]
main(tree1, tree2)


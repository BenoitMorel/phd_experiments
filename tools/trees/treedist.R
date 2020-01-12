library(ape)
library(phangorn)


get_trees <- function(treenames_file)
{
  f=file(treenames_file,open="r")
  lines=readLines(f) 
  close(f)
  return(lines)
}

main <- function(families_path, treenames_file, output_dir)
{
  families <- list.files(families_path)
  runs <- get_trees(treenames_file)
  true_tree_str <- runs[1]
  kf_vector <- vector(mode="character", length=length(runs))
  for (family in families) {

    output_path <- paste(output_dir, family, sep="/")
    print(output_path)
    prefix <- paste(families_path, family, sep="/")
    prefix <- paste(prefix, "gene_trees", sep="/")
    true_tree_path <- paste(prefix, true_tree_str, sep="/")
    true_tree <- unroot(read.tree(true_tree_path))
    index <- 1
    for (run in runs) {
      generax_path <- paste(prefix, run, sep="/")
      generax_tree <- unroot(read.tree(generax_path))
      run_kf <- mean(KF.dist(true_tree, generax_tree))
      kf_vector[index] <- run_kf
      index <- index + 1
    }
    writer <- file(output_path)
    writeLines(kf_vector, writer, sep = " ")
    close(writer)
  }
}


args = commandArgs(trailingOnly=TRUE)
families_path <- args[1]
treenames_file <- args[2]
output_dir <- args[3]


print(families_path)
main(families_path, treenames_file, output_dir)



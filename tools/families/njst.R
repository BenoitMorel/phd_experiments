library(phybase)

read_lines <- function(filename)
{
  f=file(filename,open="r") 
  lines=readLines(f) 
  close(f)
  return(lines)
}

read_first_line <- function(filename)
{
  return (read_lines(filename)[1])
}

main <- function(sptree_file, gtrees_file)
{
  sptree <- read_first_line(sptree_file)
  spname<-species.name(sptree)
  print(spname)
  nspecies<-length(spname)
  
  nodematrix<-read.tree.nodes(sptree,spname)$node
  seq<-rep(1,nspecies)
  species.structure<-matrix(0,nspecies,nspecies)
  diag(species.structure)<-1
  print(species.structure)

  genetrees = read_lines(gtrees_file)
  for(i in 1:length(genetrees))
  {
    read.tree.nodes(genetrees[i])
  }
  print(genetrees)
  res<- NJst(genetrees, spname, spname, species.structure)
  ##print(res)
}

args = commandArgs(trailingOnly=TRUE)
sptree_file <- args[1]
gtrees_file <- args[2]
main(sptree_file, gtrees_file)


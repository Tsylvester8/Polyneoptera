# helper functions

# sample the data to match tree tips and to account for uncertainty
# in tip states
GetData <- function(trees="../data/trees/posterior.trees.nex", 
                    data="../data/chrom.data/chroms.csv"){
  # get packages
  library(chromePlus)
  library(diversitree)
  
  # read trees
  trees <- read.nexus(file = trees)
  
  # read data
  chroms <- read.csv(file = data, 
                     as.is = T)
  
  
  # get exact matches
  new.dat <- as.data.frame(matrix(,1,5))
  counter <- 1
  hit.genera <- c()
  hit.family <- c()
  colnames(new.dat) <- c("order", "species", "haploid", "sex.system", "SCS")
  for(i in 1:length(trees[[1]]$tip.label)){
    current <- trees[[1]]$tip.label[i]
    if(current %in% chroms$species){
      hit <- which(chroms$species == current)
      if(length(hit) > 1) hit <- sample(hit, 1)
      new.dat[counter, ] <- chroms[hit, c(1,4,10,5,6)]
      hit.genera <- c(hit.genera, chroms[hit,3])
      hit.family <- c(hit.family, chroms[hit, 2])
      counter <- counter + 1
    }
  }
  # find genera found
  hit.genera <- unique(hit.genera)
  # find genera not found
  unhit.genera <- unique(chroms$genus[!chroms$genus %in% hit.genera])
  # split names at underscores
  tree.taxa <- strsplit(trees[[1]]$tip.label, "_")
  # making a table with genus and species epethet split out to columns
  tree.taxa <- matrix(unlist(tree.taxa), 
                      length(unlist(tree.taxa)), 
                      2, byrow=T)
  for(i in 1:length(unhit.genera)){
    if(unhit.genera[i] %in% tree.taxa[,1]){
      hit <- which(tree.taxa[,1] == unhit.genera[i])
      if(length(hit) > 1) hit <- sample(hit, 1)
      hit <- which(chroms$genus == tree.taxa[hit, 1])
      if(length(hit) > 1) hit <- sample(hit, 1)
      new.dat[counter, ] <- chroms[hit, c(1,4,10,5,6)]
      new.dat[counter, 2] <- paste(strsplit(new.dat[counter,2], "_")[[1]][1], "_sp",sep="")
      counter <- counter + 1
      hit.family <- c(hit.family, chroms[hit, 2])
    }
  }
  return(new.dat)
}
# Terrence Sylvester
# February 12 2019
# Analysis of ancestral state reconstructions

# order Orthoptera

# load required libraries
library(ape)
library(phytools)

chroms.vec <- c()

# make the plot
for(i in 1:100){
  
  # for each tree read in the probability table
  probs.input <- paste("../asr.chrom.num/Orthoptera/",
                       i, 
                       "/results/ancestorsProbs.txt", 
                       sep = "")
  
  chroms <- read.table(file = probs.input, 
                       header = T)
  
  # get the highest probable chromosome count
  max(chroms[1,-1])
  
  match(max(chroms[1,-1]), chroms[1,-1])
  
  chroms.vec[i] <- (names(chroms[-1])[match(max(chroms[1,-1]), chroms[1,-1])])
}

for( i in 1:length(unique(chroms.vec))){
  print(unique(chroms.vec)[i])
  print(sum(chroms.vec == unique(chroms.vec)[i]))
}


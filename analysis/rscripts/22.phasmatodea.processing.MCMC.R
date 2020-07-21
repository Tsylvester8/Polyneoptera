# Terrence Sylvester
# March 18 2019
# pradakshanas@gmail.com

# scaling the chromosome evolution rates of phasmatodea by tree depth

# load required libraries
library(phytools)

# read in trees
trees <- read.nexus(file = "../data/trees/phas.trees.nex") # change the data folder

# read in the chromosome rates
load(file = "../results/21.phasmatodea.MCMC.RData")

# get the tree depth for each tree and scale the rates
tree.depth <- vector()

for (i in 1:100) {
  tree.depth[i] <- max(branching.times(trees[[i]]))
  results[[i]][2:13] <- results[[i]][2:13] / tree.depth[i]
}

# discard the burnin and retain 10000 points

# user input
sample.size <- 10000
burnin <- .5

# calclating for user
start <- round(nrow(results[[1]]) * burnin)
end <- nrow(results[[1]])
l.out <- round(sample.size/length(results))

# setting up results
phas.rates <- results[[1]][0, ]

# assembling results
for(i in 1:length(results)){
  c.rows <- sample(start:end, size = l.out, replace = F)
  phas.rates <- rbind(phas.rates, results[[i]][c.rows, ])
}

# remove all but the sampled resutls
rm(list = ls()[-6])

# write the results
write.csv(phas.rates, file="../results/22.pruned.phas.rates.csv", row.names=F)


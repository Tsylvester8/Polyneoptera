# Terrence Sylvester
# March 18 2019
# pradakshanas@gmail.com

# scaling the chromosome evolution rates of phasmatodea by tree depth

# load required libraries
library(phytools)
library(ape)

# read in trees
trees <- read.nexus(file = "../data/trees/phas.trees.nex") # change the data folder

dat <- read.csv(file = "../data/chroms/phas.chroms.csv",
                as.is = T)

probs <- matrix(data = NA, nrow = 41,ncol = 2)
colnames(probs) <- c("parth", "sexual")
probs[,1] <- dat$parthenogenetic
probs[,2] <- dat$sexual

row.names(probs) <- dat$species

tree.depth <- c()

#store simmaps
simmaps <- vector(mode = "list", length = 100)
for(i in 1:100){
  tree <- keep.tip(trees[[i]], dat$species)
  
  tree.depth[i] <- max(branching.times(tree))
  
  tree$edge.length <- tree$edge.length / tree.depth[i]
  
  print(i)
  mod <- matrix(c(0,1,0,0),2,2)
  simmaps[[i]] <- make.simmap(tree = tree,
                         x = probs,
                         model =mod,
                         nsim = 100)
}

# get the number of transitions from sexual to parth
sex.to.parth <- parth.to.sex <- c()

for(i in 1:100){
  map.summary <- describe.simmap(simmaps[[i]])
  sex.to.parth[i] <- mean(map.summary$count[,3])
  parth.to.sex[i] <- mean(map.summary$count[,1])
}

# mean number of transitions from xo to xy
mean(sex.to.parth)

# mean number of transitions from xy to xo
mean(parth.to.sex)

# save reasults
save.image("../results/23.phasmatodea.stochastic.mapping.RData")

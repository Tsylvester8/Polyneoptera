# Terrence Sylvester
# June 29 2019
# pradakshanas@gmail.com

# ancestral states reconstruction for sex chromosome systems

# load helper function
source("helper.functions.R")

# load required libraries
library(phytools)

# read in tree data
trees <- read.nexus("../data/Trees/posterior.trees.nex")

# store named sex chromosome states for all trees
scs <- vector(mode = "list", length = 100)

# get orders list
# read in chromosome data
dat <- GetData(data = "../data/chrom.data/chroms.csv",
               trees = "../data/Trees/posterior.trees.nex")

#store simmaps
simmaps <- vector(mode = "list", length = 100)

for(i in 1:100){
  print(i)
  
  # read in chromosome data
  dat <- GetData(data = "../data/chrom.data/chroms.csv",
                 trees = "../data/Trees/posterior.trees.nex")
  
  # correction for Notoptera
  dat$order[dat$order %in% c("Grylloblata", "Mantophasmatodea")] <- "Notoptera"
  
  # remove taxa with no data for sex chromosome system
  dat <- dat[dat$SCS != "",]
  
  # for this ancestral states reconstruction, all XY systems (imcluding complex XY)
  # are labeled as XY. This is because transitions from XO sex chromosome system
  # to multi-XY sex chromosome system should include an intermediate 
  # sex chromosome system
  
  dat$SCS[dat$SCS %in% c("complex XY",
                         "complex XY|homomorphic",
                         "homomorphic",
                         "XY|homomorphic",
                         "XY")] <- "XY"
  
  # all other non XY sex chromosome systems are labled as XO
  dat$SCS[dat$SCS != "XY"] <- "XO"
  
  # make a named vector wich includes species name and sex chromosome system (scs)
  scs[[i]] <- setNames(dat$SCS,dat$species)
  
  # select a tree i
  tree <- trees[[i]]
  
  # remove tips that do not have data
  tree <- keep.tip(tree,dat$species)
  
  # fit ARD model to estimate the ancestral states.
  # an ARD model states unequal transisions between states
  simmaps[[i]] <- make.simmap(tree = tree,
                              x = scs[[i]],
                              model = "ARD",
                              nsim = 100)
}

xo.to.xy <- xy.to.xo <- c()

for(i in 1:100){
  map.summary <- describe.simmap(simmaps[[i]])
  xo.to.xy[i] <- mean(map.summary$count[,2])
  xy.to.xo[i] <- mean(map.summary$count[,3])
}

# mean number of transitions from xo to xy
mean(xo.to.xy)

# mean number of transitions from xy to xo
mean(xy.to.xo)

# save reasults
save.image("../results/sex.chrom.simmap.RData")

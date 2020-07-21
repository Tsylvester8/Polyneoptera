# Terrence Sylvester
# 1 April 2019

# MCMC for misbehaved trees in the order Phasmatodea

# this function allows for matching and random choice when multiple records
source("02.helper.functions.R")

# load the latest version of chromePlus
library(devtools)
install_github('coleoguy/chromePlus')
library(chromePlus)

# load diversitree for mcmc and mkn functions
library(diversitree)

# read in data to get order names
dat <- GetData(trees = "../data/trees/posterior.trees.nex", 
               data = "../data/chrom.data/chroms.csv")

# misbehaved trees
trees.to.correct <- c(27,28,39,66,81)

# to store the results
results.trees <- vector("list", length = length(trees.to.correct))

names(results.trees) <- trees.to.correct


# set values for mcmc
iter <- 1000

# read trees that will be used in analysis
trees <- read.nexus("../data/trees/posterior.trees.nex")

# counter
i <- 1

for(j in trees.to.correct[1]){
  dat <- GetData(trees = "../data/trees/posterior.trees.nex", 
                 data = "../data/chrom.data/chroms.csv")
  # drop any tips from tree without data
  tree <- drop.tip(phy = trees[[j]], tip = dat$species[is.na(dat$haploid)])
  # drop all but focus order from tree
  tree <- keep.tip(phy = tree, tip = dat$species[dat$order == "Phasmatodea"])
  tree.depth <- max(branching.times(tree))
  # scale tree to unit length
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  # create reduced dataset
  dat.t <- dat[dat$order == "Phasmatodea", ]
  # track progress
  print(paste("working on Phasmatodea tree ", j))
  # make the data format for chromePlus
  rng <- c(range(dat.t$haploid, na.rm = T)[1] - 2,
           range(dat.t$haploid, na.rm = T)[2] + 2)
  chrom.mat <- datatoMatrix(x = dat.t[!is.na(dat.t$haploid), c(2, 3)], 
                            range = rng,
                            hyper = F)
  # get the basic likelihood function 
  lik <- make.mkn(tree = tree, 
                  states = chrom.mat, 
                  k = ncol(chrom.mat),
                  strict = F,
                  control = list(method = "ode"))
  # constrain to a biologically realistic model of chrom evolution
  con.lik <- constrainMkn(data = chrom.mat, 
                          lik = lik, 
                          polyploidy = F, 
                          hyper = F,
                          constrain=list(drop.demi=T))
  # prior with little info
  # priors <- make.prior.exponential(.5)
  # sampled from an mcmc
  results.trees[[i]] <- mcmc(con.lik, w = 1,
                             x.init = runif(min = 0, max = 1, n = 3), 
                             nsteps = iter, upper = 100)
  # scale units to millions of years
  
  results.trees[[i]][,2:4] <- results.trees[[i]][,2:4] / tree.depth
  i <- i + 1
}

# plot the results and see whether these runs are behaving 
# normally or whether they need to be re run. Repeat this step until
# you get a good MCMC run for these oddly behaved trees.

for (k in 1:length(results.trees)) {
  plot(results.trees[[k]]$p,
       type = "l",
       main = paste("Tree", names(results.trees)[k]))
}

# save the results
save.image("../results/07.phasmatodea.correction.RData")


##
# clear the environment tab
rm(list = ls())

# read in the order analysis
load("../results/05b.order.analysis-4-5.RData")

# keep all but the rates
rm(list = ls()[c(-12,-18)])

# load corrected phasmatodea runs
load("../results/07.phasmatodea.correction.RData")

# replace the misbehaved runs of phasmatodea with the 
# corrected runs

for(i in trees.to.correct){
  results.orders$Phasmatodea[[i]] <- results.trees[[which(names(results.trees) == i)]]
}

# remove unwanted objects from the environmental panel
rm(list = ls()[-11])

# save the results of running this analysis
save.image("../results/07.order.rates.4-5.corrected.RData")



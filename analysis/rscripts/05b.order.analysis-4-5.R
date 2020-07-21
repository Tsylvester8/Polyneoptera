# Heath Blackmon
# 24 August 2018

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
# make table of data points
ord.tab <- as.data.frame(table(dat$order))

# pick those orders with more than 20 taxa
orders <- as.character(ord.tab$Var1)[ord.tab$Freq>20]
results.trees <- vector("list", length = 100)
results.orders <- vector("list", length = length(orders))
names(results.orders) <- orders

# set values for mcmc
iter <- 1000

# read trees that will be used in analysis
trees <- read.nexus("../data/trees/posterior.trees.nex")


for(i in 4:5){
  for(j in 1:100){
    dat <- GetData(trees = "../data/trees/posterior.trees.nex", 
                   data = "../data/chrom.data/chroms.csv")
    # drop any tips from tree without data
    tree <- drop.tip(phy = trees[[j]], tip = dat$species[is.na(dat$haploid)])
    # drop all but focus order from tree
    tree <- drop.tip(phy = tree, tip = dat$species[dat$order != orders[i]])
    tree.depth <- max(branching.times(tree))
    
    # scale tree to unit length
    tree$edge.length <- tree$edge.length/max(branching.times(tree))
    # create reduced dataset
    dat.t <- dat[dat$order == orders[i], ]
    # track progress
    print(paste("working on", orders[i]))
    # make the data format for chromePlus
    rng <- c(range(dat.t$haploid, na.rm = T)[1] - 2,
             range(dat.t$haploid, na.rm = T)[2] + 2)
    chrom.mat <- datatoMatrix(x = dat.t[!is.na(dat.t$haploid), c(2, 3)], 
                              range = rng,
                              hyper = F)
    tree.order <- c()
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
    results.trees[[j]] <- mcmc(con.lik, w = 1,
                               x.init = runif(min = 0, max = 1, n = 3), 
                               nsteps = iter, upper = 100)
    # scale units to millions of years
    results.trees[[j]][,2:4] <- results.trees[[j]][,2:4] / tree.depth
  }
  results.orders[[i]] <- results.trees
}

# save the results of running this analysis
save.image("../results/05b.order.rates4-5.RData")


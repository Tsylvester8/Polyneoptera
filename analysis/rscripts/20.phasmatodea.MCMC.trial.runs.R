# Heath Blackmon
# coleoguy@gmail.com
# 23 Jan 2018

# load libraries
library(diversitree)
library(chromePlus)

# read trees
trees <- read.nexus("../data/trees/phas.trees.nex")
# read data
chroms <- read.csv("../data/chroms/phas.chroms.csv", as.is = T)


# Droptips from the tree that are not in the trait dataset 
trees <- lapply(
  trees,
  drop.tip,
  tip = c(
    "Acanthoxyla_geisovii",
    "Acanthoxyla_huttoni",
    "Acanthoxyla_prasina",
    "Acanthoxyla_suteri"
  )
)
class(trees) <- "multiPhylo"

# this script will randomly sample a tree from the posterior distribution
# and will conduct and MCMC analysis 

# this MCMC process will allow us to get better starting values
results <- list()
for (i in 1:1000) {
  tree <- trees[[sample(1:100, 1)]]
  obs.chroms <- vector()
  for (j in 1:nrow(chroms)) {
    if (chroms$fem2N[j] != "") {
      pos.chroms <- unlist(strsplit(x = chroms$fem2N[j],
                                    split = ","))
      obs.chroms[j] <- sample(pos.chroms, 1)
    }
  }
  # at this point we have a good tree
  # a single chrom value for each species
  # we are ready to start trying to anlayze the data
  new.dat <- data.frame(chroms$species,
                        as.numeric(obs.chroms),
                        chroms$sexual)
  chrom.mat <- datatoMatrix(x = new.dat,
                            range = range(new.dat[, 2]),
                            hyper = TRUE)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  lik <- make.musse(
    tree = tree,
    states = chrom.mat,
    k = ncol(chrom.mat),
    strict = FALSE,
    control = list(method = "ode")
  )
  length(argnames(lik))
  con.lik <- constrainMuSSE(
    data = chrom.mat,
    lik = lik,
    polyploidy = F,
    constrain = list(drop.demi = T)
  )
  argnames(con.lik)
  x.init <- startVals(   
  # by doing a couple of MLEs you can 
  # come up with a better StartVals function
    n = 10,              
    min = 5,
    max = 10,
    dist = "unif"
  )
  
  mle.result <- find.mle(con.lik,
                         x.init = x.init)
  mcmc.result <- mcmc(con.lik,
                      x.init = x.init,
                      w = 1,
                      nsteps = 1000)
}

# remove unwanted results and save data
save.image("../results/20.trial.RData")
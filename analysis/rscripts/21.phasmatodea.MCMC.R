# Heath Blackmon
# coleoguy@gmail.com
# 23 Jan 2018

library(diversitree)
library(chromePlus)

# read trees
trees <- read.nexus("../data/trees/phas.trees.nex")

# read data
chroms <- read.csv("../data/chroms/phas.chroms.csv", as.is = T)

# read in trial runs

# we conducted preliminary analysis of our data. We are using these data 
# to set the starting point of our MCMC analysis
trial.runs <- read.csv("../data/trial.runs/trial.csv")[,-1]

# remove the burnin from trial runs
trial.runs <- trial.runs[c(501:1000,1501:2000), ]

results <- list()

for (i in 1:100) {
  tree <- trees[[i]]
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
                            range = c(range(new.dat[, 2])[1]-4, range(new.dat[, 2])[2]+4) ,
                            hyper = TRUE)
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  lik <- make.musse(
    tree = tree,
    states = chrom.mat,
    k = ncol(chrom.mat),
    strict = FALSE,
    control = list(method = "ode")
  )
  con.lik <- constrainMuSSE(
    data = chrom.mat,
    lik = lik,
    polyploidy = F,
    constrain = list(drop.demi = F)
  )
  x.init <- apply(X = trial.runs,
                   MARGIN = 2,
                   FUN = sample,
                   size = 1)[1:12]
   results[[i]] <- diversitree::mcmc(con.lik,
                      x.init = x.init,
                      w = apply(X=trial.runs, MARGIN=2, sd)[1:12],
                      nsteps = 5000,
                      lower = 0,
                      upper = c(rep(40,9), .00001, 40, 40),
                      prior= make.prior.exponential(r=.6))
}

# remove unwanted results and save data
save.image("../results/21.phasmatodea.MCMC.RData")
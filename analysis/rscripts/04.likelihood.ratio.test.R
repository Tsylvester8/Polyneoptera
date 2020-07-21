# Heath Blackmon
# 1 May 2019

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

# read trees that will be used in analysis
trees <- read.nexus("../data/trees/posterior.trees.nex")

anova.res <- vector(mode = "list", length = 5)
names(anova.res) <- orders

# read post burnin to get mean rates
load("../results/08.post.burnin.RData")

# remove Blattodea from post burnin since this script do not look at Blattodea
# as a whole

post.burnin <- post.burnin[-1]

for(i in 1:5){
  for(j in 1:100){
    dat <- GetData(trees = "../data/trees/posterior.trees.nex", 
                   data = "../data/chrom.data/chroms.csv")
    # drop any tips from tree without data
    tree <- drop.tip(phy = trees[[j]], tip = dat$species[is.na(dat$haploid)])
    # drop all but focus order from tree
    tree <- drop.tip(phy = tree, tip = dat$species[dat$order != orders[i]])
    tree.depth <- max(branching.times(tree))
    
    #scale tree to unit length
    tree$edge.length <- tree$edge.length/max(branching.times(tree))
    # create reduced dataset
    dat.t <- dat[dat$order == orders[i], ]
    # track progress
    print(paste("working on", orders[i]))
    print(paste("working on tree ", j))
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
    con.likwp <- constrainMkn(data = chrom.mat, 
                              lik = lik, 
                              polyploidy = F, 
                              hyper = F,
                              constrain=list(drop.demi=T))
    con.likwop <- constrainMkn(data = chrom.mat, 
                               lik = lik, 
                               polyploidy = F, 
                               hyper = F,
                               constrain=list(drop.demi=T, drop.poly=T))
    argnames(con.likwp)    
    argnames(con.likwop)    
    
    fitw <- find.mle(con.likwp, x.init=c(mean(post.burnin[[i]]$asc1*tree.depth),
                                         mean(post.burnin[[i]]$desc1*tree.depth),
                                         mean(post.burnin[[i]]$pol1)*tree.depth), method = "subplex", lower = 0, upper = 145)
    
 
        
    fitwo <- find.mle(con.likwop, x.init=c(mean(post.burnin[[i]]$asc1*tree.depth),
                                           mean(post.burnin[[i]]$desc1*tree.depth)), method = "subplex", lower = 0, upper = 145)    
    
    anova.res[[i]][[j]] <- anova(fitw,fitwo)
  }
}

# remove all but the resutls
rm(list = ls()[-1])

save.image(file = "../results/04.likelihood.ratio.test.RData")

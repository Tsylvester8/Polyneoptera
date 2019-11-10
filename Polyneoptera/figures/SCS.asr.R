# Terrence Sylvester
# April 29 2019
# pradakshanas@gmail.com

# ancestral states reconstruction for sex chromosome systems

# load helper function
source("../rscripts/helper.functions.R")

# load required libraries
library(ape)
library(beeswarm)
library(phytools)

# read in tree data
trees <- read.nexus("../data/Trees/posterior.trees.nex")

# store the root state
rootState <- c()

# store named sex chromosome states for all trees
scs <- vector(mode = "list", length = 100)

# store the asr for all trees
fitARD <- vector(mode = "list", length = 100)

# get orders list
# read in chromosome data
dat <- GetData(data = "../data/chrom.data/chroms.csv",
               trees = "../data/Trees/posterior.trees.nex")

# correction for Notoptera
dat$order[dat$order %in% c("Grylloblata", "Mantophasmatodea")] <- "Notoptera"

orders <- unique(dat$order)
ord.roots <- matrix(data = NA, nrow = 9, ncol = 101)
colnames(ord.roots) <- c("order", 1:100)
ord.roots[,1] <- orders

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
  fitARD[[i]] <- ace(x = scs[[i]], phy = tree, type = "discrete", model = "ARD")
  
  # get the probability of the root state
  rootState <- rbind(rootState, fitARD[[i]]$lik.anc[1,])
  
  for (j in 1:length(orders)) {
    ord.roots[j,(i+1)] <- (getMRCA(phy = tree, tip = dat$species[dat$order == orders[j]]) - Ntip(tree))
  }
}

# pick a random tree for plotting
randTree <- 1

tree <- trees[[randTree]]

tree <- keep.tip(phy = tree, tip = names(scs[[randTree]]))

# make a colour vector for two sex chromosome systems
cols <- setNames(c("red","blue"),levels(as.factor(scs[[randTree]])))

# plot the tree
par(mar = c(rep(0.1,4)), mai = rep(0.2,4))
plot(tree, show.tip.label = F, type = "fan")

# mark nodes
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitARD[[randTree]]$lik.anc,piecol=cols,cex=0.23)

# mark tip states
tiplabels(pie=to.matrix(scs[[randTree]][tree$tip.label],levels(as.factor(scs[[randTree]]))),
          piecol=cols,cex=0.2)

# add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

# get the propabilities of the roots of the orders
ord.probs <- matrix(data = NA, nrow = 9, ncol = 101)
colnames(ord.probs) <- c("order", 1:100)
ord.probs[,1] <- orders

for(i in 1:100){
  for (j in 1:9) {
    ord.probs[j,i+1] <- fitARD[[i]]$lik.anc[as.numeric(ord.roots[j,i+1]),1]
  }
}

for (i in 1:9) {
  print(paste(ord.probs[i,1], 100* round(mean(as.numeric(ord.probs[i,-1])),4)))
}

# distribution of rates
rates <- c()
for (i in 1:100) {
  rates <- rbind(rates, fitARD[[i]]$rates)
}

colnames(rates) <- c("Rate 1", "Rate 2")

rates <- as.data.frame(rates)

dev.off()
beeswarm(rates,
         ylab = "Rate estimate",
         xlab = "Rate",
         pch = 16,
         cex = 1,
         cex.lab = 1.5)

boxplot(rates, add = T, col = rgb(0,0,1,.1))

mean(rates$`Rate 1`)
mean(rates$`Rate 2`)

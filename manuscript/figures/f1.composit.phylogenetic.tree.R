# Terrence Sylvester
# June 21 2020
# pradakshanas@gmail.com

# load libraries
library(phytools)
library(viridis)
library(geiger)

# load helper function
source("../../analysis/rscripts/02.helper.functions.R")

# load ancestral state reconstructions for sex chromosome system
load("../../analysis/results/12.sex.chrom.asr.RData")

# read in tree data
trees <- read.nexus("../../analysis/data/Trees/posterior.trees.nex")

# get a single tree
tree <- trees[[1]]

# read in chromosome data
dat <- GetData(data = "../../analysis/data/chrom.data/chroms.csv",
               trees = "../../analysis/data/Trees/posterior.trees.nex")

# correction for Notoptera
dat$order[dat$order %in% c("Grylloblata", "Mantophasmatodea")] <- "Notoptera"

# chage Xy and complex systems to XY
dat$SCS[dat$SCS %in% c("complex XY",
                       "complex XO",
                       "complex XY|homomorphic",
                       "homomorphic",
                       "XY|homomorphic",
                       "XY")] <- "XY"

# put parthenogenetic reproductive mode in scs data
dat$SCS[dat$sex.system == "parthenogenetic"] <- "parthenogenetic"

# lable taxa whith no data for scs as unknown
dat$SCS[dat$SCS == ""] <- "unknown"

# sort the data table so that the order of the species is same as that of the 
# tree
dat.ordered <- as.data.frame(matrix(data = NA,
                                    nrow = nrow(dat),
                                    ncol = ncol(dat)))

colnames(dat.ordered) <- colnames(dat)

for(i in 1:Ntip(tree)){
  dat.ordered[i,] <- dat[dat$species == tree$tip.label[i],]
}

# get the root of all the orders
ord.roots <- matrix(data = NA, nrow = 10, ncol = 2)
colnames(ord.roots) <- c("order", "root")
ord.roots[,1] <- c(orders,"origin")

for (j in 1:length(orders)) {
  ord.roots[j,2] <- (getMRCA(phy = trees[[1]], 
                             tip = dat$species[dat$order == orders[j]]) - Ntip(trees[[1]]))
}

# get the root of the tree
ord.roots[10,2] <- 1

# get the ancestral states of the orders
asr.orders <- matrix(data = NA, nrow = 10,
                     ncol = 2)
colnames(asr.orders) <- c("XO", "XY")

for (i in 1:9) {
  asr.orders[i,1] <- round(mean(as.numeric(ord.probs[i,-1])),4) * 100
  asr.orders[i,2] <- 100 - asr.orders[i,1]
}

# get the ancestral state at the root of the tree
asr.orders[10,1] <- round(mean(rootState[,1]),4) * 100
asr.orders[10,2] <- 100 - asr.orders[10,1]


# plot the tree
plotTree.wBars(tree = tree,
               x = setNames(dat.ordered$haploid, dat.ordered$species),
               type = "fan",
               scale = 4,
               method = "plotTree",
               width = 5,
               col = rainbow(9, alpha = 0.5)[as.factor(dat.ordered$order)],
               offset = 1,
               border = NA,
               lwd = 1)

# mark nodes
#set node colours
cols <- setNames(c("#4daf4a","#984ea3"), levels(as.factor(c("XO", "XY"))))

nodelabels(node=as.numeric(ord.roots[,2]) + Ntip(trees[[1]]),
           pie= asr.orders,
           piecol=cols,
           cex=.3)

# mark tips
tiplabels(pch = 16,
          col = c("#e41a1c",
                  "#377eb8",
                  "#4daf4a",
                  "#984ea3")[as.factor(dat.ordered$SCS)],
                  cex = .7)

# legend
text(x = -620, y = 525, labels = "Orders",
     pos = 4)

points(x = rep(-600, 9),
       y = seq(from = 500, to = 300, length.out = 9),
       pch = 22,
       bg =rainbow(9, alpha = .5),
       col = "black",
       cex = 1.5)

text(x = rep(-600, 9),
     y = seq(from = 500, to = 300, length.out = 9),
     labels = levels(as.factor(dat.ordered$order)),
     pos = 4,
     cex = .8)

text(x = -620, y = 250, labels = "SCS",
     pos = 4)

points(x = rep(-600, 4),
       y = c(225,200,175,150),
       pch = 22,
       bg = c("#e41a1c","#377eb8","#4daf4a","#984ea3"),
       col = "black",
       cex = 1.5)

text(x = rep(-600, 4),
     y = c(225,200,175,150),
     labels = levels(as.factor(dat.ordered$SCS)),
     pos = 4,
     cex = .8)


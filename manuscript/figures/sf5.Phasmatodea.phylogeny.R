# Terrence Sylvester
# 8 April 2019
# pradakshanas@gmail.com

# phasmatodea tree and traits

# load required packages
library(phytools)
library(viridis)
library(evobiR)

# read in data
dat <- read.csv("../data/chroms/phas.chroms.csv", as.is = T)

# read in trees
trees <- read.nexus("../data/trees/phas.trees.nex")

# Droptips from the phylogeny as these are not included in our trait
# dataset
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

# vector of sex chromosome system
SCS <- dat$XY
names(SCS) <- dat$species

# vector of chromosome number
chroms <- dat$fem2N
names(chroms) <- dat$species

# vector of mode of reproduction
parths <- dat$parthenogenetic
names(parths) <- dat$species

# plot sex chromosome system
ShowTree(tree = trees[[1]], 
         tip.vals = SCS, 
         pch = 16, 
         cols = viridis(3), 
         tip.cex = 1)

# plot a tree with chromosome numbers
cn <- vector() # store mean chromosome number

for(i in 1:length(chroms)){
  cn[i] <- mean(as.numeric(unlist(strsplit(chroms[i], 
                                           split = ","))))
}

names(cn) <- names(chroms)
tip.col <- viridis(101)[100*parths + 1]

# set the plot window
par(mar = c(0, 0, 0, 5))

# plot the tree
plot(trees[[1]], show.tip.label = F,
     type = "fan", edge.width = 3)

# plot tip lables

# colour them according to their reproductive mode
tiplabels(pch = 16, col = tip.col, cex = 2)

# mean chromosome number for each tip
tiplabels(text = round(cn,digits = 2), 
          frame = "none", cex = 1, offset = 6)

# plot the contMap of traits
phascont <-  contMap(tree = trees[[1]],
                     x = parths,
                     res = 100,
                     lwd = 4,
                     ftype = "off",
                     legend = T,
                     outline = F,
                     type = "fan",plot = F)

# change colours of the contMap
phascont$cols <- viridis(1001)
names(phascont$cols) <- 0:1000

plot.contMap(phascont,
             res = 100,
             lwd = 3,
             ftype = "off",
             outline = F,
             type = "fan")

tree

plotSimmap(tree = trees[[1]],colors = tip.col)


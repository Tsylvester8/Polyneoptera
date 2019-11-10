# figure showing the phylogenetic relationships and the haploid chromosome
# count of the taxa

# load the required libraries
library(phytools)
library(viridis)
library(diversitree)

source(file = "../rscripts/helper.functions.R")

dat <-  GetData(trees = "../data/Trees/posterior.trees.nex",
                data = "../data/chrom.data/chroms.csv")

trees <- read.nexus(file = "../data/Trees/posterior.trees.nex")

# The orders Grylloblattodea and Mantophasmatodea in
# our dataset are concidered as suborders.
# These two sub orders fall under one order called
# Notoptera

dat$order[dat$order %in% c("Grylloblata", "Mantophasmatodea")] <- "Notoptera"

# make a named vector containing the phenotypic trait
# tobe ploted into the tree
# in this case it is the haploid chromosome number

haploid.chroms <- vector()

for(i in 1:length(trees[[1]]$tip.label)){
  haploid.chroms[i] <- dat$haploid[trees[[1]]$tip.label[i] == dat$species]
}

names(haploid.chroms) <- trees[[1]]$tip.label

# plot the tree 
# log transform the haploid chromosome numbers
chrom.map <- contMap(tree = trees[[1]],
        x = log(haploid.chroms),
        lwd = 2,
        outline = F,
        type = "fan",
        res = 1000,
        ftype = "off",
        fsize = c(0.5,0.7),
        legend = 100, 
        plot = F)

# changes the colour scheme of the tree using viridis
chrom.map$cols <- viridis(n = 1001)
names(chrom.map$cols) <- 0:1000

# adgist the limits of the contMap so that it represents the real values
# not the log transformed values
chrom.map$lims <- c(min(haploid.chroms), max(haploid.chroms))

# plot the haploid chromosome numbers as tiplabels
# tiplabels(text = haploid.chroms,
#           frame = "none",
#           cex = 0.3,
#           offset = 10)

plot.contMap(chrom.map, 
             type = "fan", 
             ftype = "off",
             fsize = c(0.5,.8),
             outline = F,
             lwd = 3,
             legend = F,
             res = 1000,
             mar = rep(.5,4))

add.color.bar(200,
              chrom.map$cols,
              title="Chromosome number",
              lims=chrom.map$lims,
              digits=3,
              prompt=FALSE,
              x=-550,
              y=-350,
              lwd=4,
              fsize=1,
              subtitle="")

# mark the clades

for(i in 1:length(unique(dat$order))){
  if(length(dat$species[dat$order == unique(dat$order)[i]]) > 1 &
     unique(dat$order)[i] != "Isoptera" & unique(dat$order)[i] != "Blattodea"){
    arc.cladelabels(text = unique(dat$order)[i],
                    node = findMRCA(tree = trees[[1]],
                                    c(dat$species[dat$order == unique(dat$order)[i]])),
                    mark.node = F,
                    cex = 0.8,
                    orientation = "horizontal",
                    ln.offset = 1.02,
                    lwd = 3,
                    lab.offset =1.03)
  }
  if(length(dat$species[dat$order == unique(dat$order)[i]]) == 1){
    arc.cladelabels(text = unique(dat$order)[i],
                    node = which(trees[[1]]$tip.label == dat$species[dat$order == unique(dat$order)[i]]),
                    mark.node = F,
                    cex = 0.8,
                    orientation = "horizontal")
  } 
  
  # for Blattodea
  arc.cladelabels(text = NULL,
                  node = findMRCA(tree = trees[[1]],
                                  c(dat$species[dat$order == "Blattodea"])),
                  mark.node = F,
                  cex = 0.8,
                  orientation = "horizontal",
                  ln.offset = 1.02,
                  lwd = 3) 
  
  # for Isoptera
  arc.cladelabels(text = NULL,
                  node = findMRCA(tree = trees[[1]],
                                  c(dat$species[dat$order == "Isoptera"])),
                  mark.node = F,
                  cex = 0.8,
                  orientation = "horizontal",
                  ln.offset = 1.05,
                  lwd = 3,
                  lab.offset =1.06) 
}

# Terrence Sulvester
# pradakshanas@gmail.com
# 14 May 2020
# genome size vs tip rates of chromosome number change

# load libraries
library(ape)
library(phylolm)

# use the following helper function
ReorderData <- function(tree, data, taxa.names="row names"){
  new.data <- data
  if(is.vector(data)){
    for(i in 1:length(tree$tip.label)){
      new.data[i] <- data[which(names(data) == tree$tip.label[i])]
      names(new.data)[i] <- names(data)[which(names(data) == tree$tip.label[i])]
    }
  }
  if(is.data.frame(data) || is.matrix(data)){
    if(taxa.names == "row names"){
      row.names(new.data) <- 1:length(tree$tip.label)
      for(i in 1:length(tree$tip.label)){
        new.data[i,] <- data[which(row.names(data) == tree$tip.label[i]),]
        row.names(new.data)[i] <- row.names(data)[which(row.names(data) == tree$tip.label[i])]
      }
    }
    if(is.numeric(taxa.names)){
      for(i in 1:length(tree$tip.label)){
        new.data[i,] <- data[which(data[,taxa.names] == tree$tip.label[i]),]
      }
    }
  }
  return(new.data)
}



# read in genome size data
gsz <- read.csv("../data/genome.size/polyneopteragsz.csv", as.is = T)

# read in order names
orders <- dir(path = "../asr.chrom.num/")[c(2,7,8,9,10)]

# make a list to store reslus for all orders
anc.probs <- vector(mode = "list", length = 5)
names(anc.probs) <- orders

for(i in 1:5){
  print(orders[i])
  for(j in 1:100){
    print(j)
    # read in ancestral states
    anc.state <- read.table(paste("../asr.chrom.num/",
                                  orders[i],
                                  "/",
                                  j,
                                  "/results/ancestorsProbs.txt",
                                  sep = ""),
                            as.is = T,
                            header = T)
    # read in trees
    tree <- read.tree(paste("../asr.chrom.num/",
                            orders[i],
                            "/",
                            j,
                            "/results/allNodes.tree",
                            sep = ""))

    # get the ancestral nodes - we look at the number of tips in the tree. and then
    # we look at the edge table and see wich edges ends in 1 : Ntips. (species are
    # labled from 1 to the given number of species and nodes are labled starting from
    # 1 + total number of species).
    nTips <- length(tree$tip.label)
    nodeNum <- tree$edge[,1][which(tree$edge[,2] %in% c(1:nTips))] - nTips

    # make a data table to hold species names chromosome number at the tip,
    # chromosome number of the ancestor and the rate of change
    chromChange <- as.data.frame(matrix(data = NA,
                                        nrow = nTips,
                                        ncol = 4))
    colnames(chromChange) <- c("species",
                               "tipValue",
                               "ancestralValue",
                               "rateOfChange")

    chromChange$species <- tree$tip.label

    # get the chromosome numbers
    chromValues <- gsub(pattern = "X",
                        replacement = "",
                        x = colnames(anc.state[,-1]))

    for(k in 1:nTips){
      hit.anc <- which(anc.state$NODE %in% paste("N", nodeNum, sep = "")[k])
      hit.tip <- which(anc.state$NODE %in% tree$tip.label[k])
      chromChange$ancestralValue[k] <-  as.numeric(chromValues[anc.state[hit.anc,-1] == max(anc.state[hit.anc,-1])])
      chromChange$tipValue[k] <- as.numeric(chromValues[anc.state[hit.tip,-1] == max(anc.state[hit.tip,-1])])
      brLength <- branching.times(tree)[which(names(branching.times(tree)) == paste("N", nodeNum, sep = "")[k])]
      chromChange$rateOfChange[k] <- (chromChange$tipValue[k] - chromChange$ancestralValue[k]) / brLength
    }
    anc.probs[[i]][[j]] <- chromChange
  }
}

# sort the order of the species present
for(i in 1:5){
  for(j in 1:100){
    anc.probs[[i]][[j]] <-  anc.probs[[i]][[j]][order(anc.probs[[i]][[j]]$species),]
  }
}

for(i in 1:5){
  for(j in 1:100){
    anc.probs[[i]][[1]] <- cbind(anc.probs[[i]][[1]], anc.probs[[i]][[j]]$rateOfChange)
  }
}

# now rowbind these
anc.probs.all <- c()

for(i in 1:5){
  anc.probs.all <- rbind(anc.probs.all, anc.probs[[i]][[1]])
}

# remove unwanted columns
anc.probs.all <- as.data.frame(anc.probs.all[,-c(2,3,4)])
colnames(anc.probs.all) <- c("species",
                             paste("tree", 1:100, sep = ""))

# species level matchs
species.matchs <- anc.probs.all$species[anc.probs.all$species %in% gsz$Species]

# genera level matches
gsz.gn <- gsz[!(gsz$Species %in% species.matchs),]

# rename species name so that they reflect genera
gsz.gn$Species <- paste(gsz.gn$Genus, "_sp", sep = "")
genera.matchs <- anc.probs.all$species[anc.probs.all$species %in% gsz.gn$Species]
all.matches <- c(species.matchs, genera.matchs)

# make a new data table
dat <- as.data.frame(matrix(data = NA,
                            nrow = length(all.matches),
                            ncol = 3))
colnames(dat) <- c("species",
                   "meanRate",
                   "gsz")
dat$species <- all.matches

for (i in 1:nrow(dat)) {
  dat$meanRate[i] <- mean(as.numeric(anc.probs.all[anc.probs.all$species == dat$species[i],-1]))
  x <- (as.numeric(gsz$Mbp[gsz$Species == dat$species[i]]))
  if(length(x) == 0) x <- as.numeric(gsz.gn$Mbp[gsz.gn$Species == dat$species[i]])
  dat$gsz[i] <- mean(x)
}

# regression (directional)
fit <- lm(dat$meanRate ~ dat$gsz)
fit.sum <- summary(fit)

# regression (absolute)
fit.abs <- lm(abs(dat$meanRate) ~ dat$gsz)
fit.abs.sum <- summary(fit.abs)

# phylogenetically corrected regression
phy <- read.nexus("../data/Trees/posterior.trees.nex")[[sample(1:100, 1)]]
phy <- keep.tip(phy =phy, tip = dat$species)

# reorder data
dat <- ReorderData(tree = phy, data = dat, taxa.names = 1)
rownames(dat) <- dat$species

# phylogenetic regression - directional
phylo.fit <- phylolm(formula =meanRate ~ gsz, data = dat,phy = phy)
phylo.fit.sum <- summary(phylo.fit)

# phylogenetic regression - absolute
phylo.fit.abs <- phylolm(formula =abs(meanRate) ~ gsz,data = dat,phy = phy)
phylo.fit.abs.sum <- summary(phylo.fit.abs)

# keep only results that are needed
rm(list = ls()[-c(8,11,12,27,28)])

# save results
save.image("../results/18.gsz.vs.tip.rate.RData")

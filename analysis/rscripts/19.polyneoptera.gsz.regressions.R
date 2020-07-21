# Terrence Sylvester
# 20th February 2020
# pradakshanas@gmail.com

# load libraries
library(phytools)
library(phylolm)
library(evobiR)

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
gsize <- read.csv("../data/genome.size/polyneopteragsz.csv", as.is = T)

# read in chromosome number data
chroms <- read.csv("../data/chrom.data/chroms.csv", as.is = T)

# read in trees
trees <- read.nexus("../data/Trees/posterior.trees.nex")

# fill spaces in the species names to match with that of the phylogeny
gsize$Species <-  gsub(pattern = " ",
                       replacement = "_",
                       x = gsize$Species)

# make sure that the chromosome number are numeric
gsize$Chrom.num <- as.numeric(gsize$Chrom.num)

# get the hapliod chromosome number 
gsize$Chrom.num <- gsize$Chrom.num / 2


# we find that in that dataset downloaded from the genome size database there 
# two some species whose order name is the old version of the current name. 
# fix this

gsize$Order[gsize$Order == "Blattaria"] <- "Blattodea"
gsize$Order[gsize$Order == "Phasmida"] <- "Phasmatodea"

# when there are multiple records for genome size for a given species we will
# use the mean genome size

# get the speceis for which there are more than one record 
gsize.dups <- gsize$Species[(duplicated(gsize$Species) | duplicated(gsize$Species, fromLast = T))]

# get the unique names of these duplicated species
gsize.dups.unique <- unique(gsize.dups)

# make a matrix to store the mean genome size of these duplicated records
gsize.dups.mean <- as.data.frame(matrix(data = NA, 
                          nrow = length(gsize.dups.unique),
                          ncol = ncol(gsize)))

# name the columns
colnames(gsize.dups.mean) <- colnames(gsize)

# fill in spepcies names
gsize.dups.mean$Species <- gsize.dups.unique

# fill in other information
for(i in 1:length(gsize.dups.unique)){
  gsize.dups.mean$Order[i] <- unique(gsize$Order[gsize$Species == gsize.dups.unique[i]])
  gsize.dups.mean$C.value[i] <- sample(gsize$C.value[gsize$Species == gsize.dups.unique[i]],1)
  gsize.dups.mean$Mbp[i] <- sample(gsize$Mbp[gsize$Species == gsize.dups.unique[i]],1)
  hit <- (unique(gsize[gsize$Species == gsize.dups.unique[i], -c(8,9)]))
  gsize.dups.mean[i,-c(8,9)] <- hit[1,]
}

# combine this matrix with the uniue records in the genome size dataset
gsize <- rbind(gsize.dups.mean, gsize[!(gsize$Species %in% gsize.dups),])

# get the overlapping species names between the two datasets
sp.names <- chroms$species[(chroms$species %in% gsize$Species)]

# find duplicated names in genome size database
gsize.sp <- gsize$Species[gsize$Species %in% sp.names]

hit <- gsize.sp[duplicated(gsize.sp)]

# remove them
sp.names <- sp.names[!(sp.names %in% hit)]

# fill in chromosome numbers
for(i in 1:length(sp.names)){
  if(!is.na(gsize$Species[gsize$Species %in% sp.names[i]])){
    x <- chroms$hap[chroms$species==sp.names[i]]
    if(length(x) == 1) gsize$Chrom.num[gsize$Species %in% sp.names[i]] <- x
    if(length(x) > 1) gsize$Chrom.num[gsize$Species %in% sp.names[i]] <- sample(x, 1)
  }
}

# remove species with no chromosome number information from our genome size data set
gsize.sp.level <- gsize[!is.na(gsize$Chrom.num),]

# upto now we have all species-level matches. Here we will get genus-level matches
# get records where we do not have chromosome number data from the genome size data set
gsize.no.chrom.data <- gsize[is.na(gsize$Chrom.num),]

# get the chromosome data for those who don't have species level match
chrom.gn.level <- chroms[!(chroms$species %in% gsize.sp.level$Species),]

# in this dataset we find some ambiguous genus names such as "genus7", "genus8"
# "unknown" and "Sp.". remove these

chrom.gn.level <- chrom.gn.level[-c(which(chrom.gn.level$genus %in% c("Genus7","Genus8","unknown","Sp."))),]

# get genus level chrom count
# uniqe genera
unique.genera <- unique(chrom.gn.level$genus)

chrom.gn.level.matches <- matrix(data = NA,
                                 nrow = length(unique.genera),
                                 ncol = ncol(chrom.gn.level))

colnames(chrom.gn.level.matches) <- colnames(chrom.gn.level)

chrom.gn.level.matches <- as.data.frame(chrom.gn.level.matches)

for (i in 1:length(unique.genera)) {
  chrom.gn.level.matches$species[i] <- paste(unique.genera[i], "_sp", sep = "")
  chrom.gn.level.matches$order[i] <- unique(chrom.gn.level$order[chrom.gn.level$genus == unique.genera[i]])
  chrom.gn.level.matches$hap[i] <- sample(chrom.gn.level$hap[chrom.gn.level$genus == unique.genera[i]], size = 1)
}

# get the genus level matches for genome size data
gsize.unique.genera <- unique(gsize.no.chrom.data$Genus)

gsize.gn.level.matches <- matrix(data = NA,
                                 nrow = length(gsize.unique.genera),
                                 ncol = ncol(gsize.no.chrom.data))

colnames(gsize.gn.level.matches) <- colnames(gsize.no.chrom.data)

gsize.gn.level.matches <- as.data.frame(gsize.gn.level.matches)

for (i in 1:length(gsize.unique.genera)) {
  gsize.gn.level.matches$Species[i] <- paste(gsize.unique.genera[i], "_sp", sep = "")
  gsize.gn.level.matches$Order[i] <- unique(gsize.no.chrom.data$Order[gsize.no.chrom.data$Genus == gsize.unique.genera[i]])
  gsize.gn.level.matches$Mbp[i] <- sample(gsize.no.chrom.data$Mbp[gsize.no.chrom.data$Genus == gsize.unique.genera[i]], size = 1)
}

# fill in chromosome numbers for genus-level matchs
for(i in 1:nrow(gsize.gn.level.matches)){
  hit <- which(chrom.gn.level.matches$species %in% gsize.gn.level.matches$Species[i])
  if(length(hit) == 1){
    gsize.gn.level.matches$Chrom.num[i] <- chrom.gn.level.matches$hap[hit]
  }
}

# Remove genus level matchs where we do not have chromosome data
gsize.gn.level.matches <- gsize.gn.level.matches[!(is.na(gsize.gn.level.matches$Chrom.num)),]

# make the final dataset
gsize <- rbind(gsize.gn.level.matches, gsize.sp.level)

# keep only columns tha are neded
gsize <- gsize[,c(4,7,9,10)]

#### linear models ####
### phylogenetic lm

# get a random single tree
tree <- trees[[sample(x = 1:100, size = 1)]]

tree <- keep.tip(phy = tree, tip = gsize$Species[gsize$Species %in% tree$tip.label])

# remove the data from the genome size dataset that are absent in the tree
gsize.phylo.lm <- gsize[gsize$Species %in% tree$tip.label,]

# linear model not corrected for phylogeny
gsz.lm.sum <- summary(lm(gsize$Chrom.num~gsize$Mbp))

# make a new dataset in which the data is sorted according to the order in the 
# tree
gsize.sorted <-  ReorderData(tree = tree,
                             data = gsize.phylo.lm,
                             taxa.names=2)

# rename the rows so that rownames are species names
# this is for phylogenetic regressions
rownames(gsize.sorted) <- gsize.sorted$Species

# linear model corrected for phylogeny
phylo.gsz.lm <- phylolm(formula = Chrom.num ~ Mbp, 
               data = gsize.sorted,
               phy = tree, model = "BM")

# get the summary statistics
phylo.gsz.lm.sum <- summary(phylo.gsz.lm)

# clear unwanted data
rm(list = ls()[-c(4,11,15,19)])

# save 
save.image("../results/19.polyneoptera.gsz.regressions.RData")

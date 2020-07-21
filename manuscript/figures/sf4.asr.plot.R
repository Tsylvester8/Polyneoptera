# load required libraries
library(ape)
library(phytools)
library(viridis)

# get the orders list
orders <- dir(path = "../mod.comp/")[c(-3,-4,-5,-6,-11)]

# rename Blattodea accordingly
orders.corrected <- orders

orders.corrected[c(1,2)] <- c("Blattodea Including Isoptera",
                              "Blattodea excluding Isoptera")

# set the plot pane parameters accordingly
par(mfrow = c(6,2),
    mar=c(3,1,1.5,6))

# make the plot
for(i in 1:length(orders)){
  
  # to hold the possible chromosomes at the root
  chroms.vec <- c()
  
  # to hold the probabilitites of the possible chromosomes at
  # the root
  probs.vec <- c()
  
  # pick a random tree
  tree.input <- paste("../mod.comp/", 
                      orders[i], "/", 
                      sample(1:100, size = 1), 
                      "/results/mlAncestors.tree", 
                      sep = "")
  
  # read in the tree
  tree <- read.newick(file = tree.input)
  
  #tree$edge.length <-  tree$edge.length / max(branching.times(tree))
  
  # Get the estimates of the chromosome number at the root
  # taken from the ChromEvol analysis for all the
  # 100 trees
  for(j in 1:100){
    
    # for each tree read in the probability table
    probs.input <- paste("../mod.comp/", 
                         orders[i],
                         "/", 
                         j, 
                         "/results/ancestorsProbs.txt", 
                         sep = "")
    
    chrom.probs <- read.table(file = probs.input, 
                              header = T)
    
    # get the chromosome numbers and there probabilitites
    
    # set the chrom.count acoordingly so that you get the
    # top x posible chromosome number for each tree and
    # their probabilities
    chrom.count <- 2
    maxChroms <- names(sort(chrom.probs[1,-1],decreasing = T))
    maxChroms <- as.numeric(sub("X", "", maxChroms))
    chroms.vec <- c(chroms.vec, maxChroms)
    maxProbs <- as.character(sort(chrom.probs[1,-1],decreasing = T))
    probs.vec <- c(probs.vec, maxProbs)
  }
  
  # plot the tree
  plot(tree, show.tip.label = F,
       cex = .4,label.offset = .5,
       edge.width = 1,
       node.depth = 1)
  
  # name the tree by its order
  mtext(text = paste(letters[i], ") ", orders.corrected[i], sep = ""),
        side = 3, line = 0,
        cex = .8,
        adj = 0)
  
  # lable the root node
  # nodelabels(node = tree$edge[1,1],
  #            pch = 16,
  #            col = "red")
  
  # get the range of the chromosome numbers to set the 
  # X axis limits
  rng <- range(chroms.vec)
  
  # plot the probabilities 
  plot(x = jitter(chroms.vec,factor = .4),
       y = probs.vec,
       xlim = rng,
       ylim = c(0,1),
       ylab = NA,
       xlab = NA,
       pch = 16,
       cex = .3,
       col= rgb(0,0,1,.5))
  
  # label the Y axis
    mtext(text = "Probability",side = 2, line = 4,
          cex = .7)
}

# label the X axis
mtext(text = "Number of chromosomes",side = 1, line = 2,
      cex = .7)



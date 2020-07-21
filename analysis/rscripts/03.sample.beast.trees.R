# T Sylvester 
# 11th August 2018
# pradakshanas@gmail.com

# load the required packages
library(ape)

# random sampling of trees from each BEAST run
trees <- read.nexus("../data/Trees/beast/run1.nexus")[sample(501:1000, 50)]
trees[51:100] <- read.nexus("../data/Trees/beast/run2.nexus")[sample(501:1000, 50)]
write.nexus(trees, file="../data/Trees/posterior.trees.nex")
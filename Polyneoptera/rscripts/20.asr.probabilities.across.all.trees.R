# Terrence Sylvester
# July 01 2019
# pradakshanas@gmail.com

# get the orders list
orders <- dir(path = "../mod.comp/")[c(-3,-4,-5,-6,-11)]

# following code will get the probabilities for chromosome numbers
# for each tree

# make the input file which locates the table with the probabilities
# for each chromosome numer

results <- vector(mode = "list", length = length(orders))

names(results) <- orders

for(i in 1:length(orders)){
  chrom.mat <- matrix(data = NA, 
                      nrow = 0, 
                      ncol = 102)
  
  chrom.mat <- as.data.frame(chrom.mat)
  
  colnames(chrom.mat) <- c("chrom.num", 1:100, "probability")
  
  for(j in 1:100){
    input <-  paste("../mod.comp/", 
                    orders[i],
                    "/", 
                    j, 
                    "/results/ancestorsProbs.txt", 
                    sep = "")
    
    chrom.probs <- read.table(file = input, 
                              header = T)
    
    if(nrow(chrom.mat) == 0){
      chrom.mat[nrow(chrom.mat)+length(colnames(chrom.probs)[-1]),] <- NA  
    }
    chrom.mat[,1] <- colnames(chrom.probs)[-1]  
    chrom.mat[,j+1] <- as.numeric(chrom.probs[1,-1])
  }
  results[[i]] <- chrom.mat
}

for(i in 1:length(results)){
  for(j in 1:nrow(results[[i]])){
    results[[i]]$probability[j] <- sum(results[[i]][j,2:101])/100
  }
}

# remove unwanted results and save
save.image(file = "../results/asr.probabilities.across.all.trees.RData")

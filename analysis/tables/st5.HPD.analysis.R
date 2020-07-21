# Terrence Sylvester
# 11 February 2019
# Comparison of chromosome rates for the order level analysis

# HPD intervals

# load required libraries
library(coda)

options(scipen=999)

# load the post burnin from the MCMC run
load(file = "../results/post.burnin.RData")

orders <- names(post.burnin)

# make a matrix to hold the values
HPD.tab <- matrix(data = NA,
                  nrow = length(orders),
                  ncol = 4)

colnames(HPD.tab) <- c("order",
                       "Fission",
                       "Fusion",
                       "Polyploidy")

HPD.tab <- as.data.frame(HPD.tab)

HPD.tab$order <- orders

digits <- 3

for (i in 1:length(orders)) {
  #Fissions
  HPD.tab$Fission[i] <-  paste(round(HPDinterval(as.mcmc(post.burnin[[i]]$asc1))[1],digits), 
                                     " - ",
                                     round(HPDinterval(as.mcmc(post.burnin[[i]]$asc1))[2],digits),
                                     sep = "")
  #Fusions
  HPD.tab$Fusion[i] <-  paste(round(HPDinterval(as.mcmc(post.burnin[[i]]$desc1))[1],digits), 
                              " - ",
                              round(HPDinterval(as.mcmc(post.burnin[[i]]$desc1))[2],digits),
                              sep = "")
  
  #Polyplody
  HPD.tab$Polyploidy[i] <-  paste(round(HPDinterval(as.mcmc(post.burnin[[i]]$pol1))[1],digits), 
                                  " - ",
                                  round(HPDinterval(as.mcmc(post.burnin[[i]]$pol1))[2],digits),
                                  sep = "")
  
}

# save the HPD table
write.csv(x = HPD.tab,
          file = "../tables/HPD.Table.csv",row.names = F)

# Terrence Sylvester
# June 30 2019
# pradakshanas@gmail.com

# read in data
load("../results/08.post.burnin.RData")

# make a matrix
mean.est <- matrix(data = NA,
                   nrow = length(post.burnin),
                   ncol = 4)

colnames(mean.est) <- c("order", "asc", "desc", "pol")

mean.est[,1] <- names(post.burnin)

# get means of rate estimates
for(i in 1:length(post.burnin)){
  mean.est[i,2] <- round(mean(post.burnin[[i]][[2]]),3)
  mean.est[i,3] <- round(mean(post.burnin[[i]][[3]]),3)
  mean.est[i,4] <- round(mean(post.burnin[[i]][[4]]),3)
}



# Terrence Sylvester
# June 30 2019
# pradakshanas@gmail.com

# make csv file of the order rates

# read in data
load("../results/post.burnin.RData")

# make an empty matrix
dat <- matrix(data = NA,
              nrow = 5,
              ncol = 0)

# fill the empty matrix with data
for(i in 1:length(post.burnin)){
mat <- cbind(rep(names(post.burnin)[i], nrow(post.burnin[[i]])),
      post.burnin[[i]][2:5])

dat <- rbind(dat, mat)
}

colnames(dat) <- c("order", "asc", "desc", "pol", "lnL")

# write results
write.csv(x = dat,
          file = "../results/order.rates.post.burnin.csv",row.names = F)

# Terrence Sylvester
# July 07 2019
# pradakshans@gmail.com

# load results from stochastic mapping (Phasmatodea)
load("../results/phasmatodea.stochastic.mapping.RData")

# load required libraries
library(coda)

# get the mean rate for sex to parth transitions
rate.sex.to.path <- c()

for(i in 1:100){
  for(j in 1:100){
    rate.sex.to.path <-c(rate.sex.to.path, simmaps[[i]][[j]]$Q[2,1] / tree.depth[i])
  }
}

# get the mean rate of transion s
mean.rate.sex.to.path <-  mean(rate.sex.to.path)

rate.sex.to.path.HPD <- HPDinterval(as.mcmc(rate.sex.to.path))

sex.to.parth.tab <- matrix(data = NA, nrow = 1, ncol = 2)

colnames(sex.to.parth.tab) <- c("Mean rate (95% Credible Interval)", 
                                "Mean number of transitions")

sex.to.parth.tab[1,1] <- paste(round(mean.rate.sex.to.path, 4),
                               " (",
                               round(rate.sex.to.path.HPD[1],4),
                               " - ",
                               round(rate.sex.to.path.HPD[2],4),
                               ")",
                               sep = "")

sex.to.parth.tab[1,2] <- round(mean(sex.to.parth),1)

write.csv(x = sex.to.parth.tab,
          file = "../tables/phasmatodea.stochastic.mapping.summary.csv", 
          row.names = F)

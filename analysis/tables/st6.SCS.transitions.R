load("../results/sex.chrom.simmap.RData")

# get the rates of transitions

rate.xo.to.xy <- rate.xy.to.xo <- c()
for(i in 1:100){
  for(j in 1:100){
    rate.xo.to.xy <- c(rate.xo.to.xy, simmaps[[i]][[j]]$Q[1,2])
    rate.xy.to.xo <- c(rate.xy.to.xo, simmaps[[i]][[j]]$Q[2,1])
  }
}

# get the 95% credible intervals 
library(coda)

xo.to.xy.HPD <- HPDinterval(as.mcmc(rate.xo.to.xy))[1,]
xy.to.xo.HPD <- HPDinterval(as.mcmc(rate.xy.to.xo))[1,]

# make the table
summary.tab <- matrix(data = NA, nrow = 2, ncol = 3)

colnames(summary.tab) <- c("Transition",
                           "Mean rate (95% credible interval)",
                           "Mean number of transitions")

summary.tab[,1] <- c("XO to XY",
                     "XY to XO")

summary.tab[1,2] <- paste(round(mean(rate.xo.to.xy), 4),
                          " (", 
                          round(xo.to.xy.HPD[1],4), 
                          " - ",
                          round(xo.to.xy.HPD[2],4),
                          ")",
                          sep = "")

summary.tab[2,2] <- paste(round(mean(rate.xy.to.xo), 4),
                          " (", 
                          round(xy.to.xo.HPD[1],4), 
                          " - ",
                          round(xy.to.xo.HPD[2],4),
                          ")",
                          sep = "")
                          
summary.tab[1,3] <- round(mean(xo.to.xy), 4)

summary.tab[2,3] <- round(mean(xy.to.xo), 4)

# save the table
write.csv(x = summary.tab,
          file = "../tables/SCS.transitions.csv",
          row.names = F)

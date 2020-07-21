# Terrence Sylvester
# 1st April 2019

# Analysis of MCMC

## read in data
# load the Blat-iso-combined rates
load("../results/05c.order.analysis-blat.iso.combined.RData")

# rename results accordingly
Blat.Iso.comb <- results

# remove all but the rates
rm(list = ls()[-1])

# look at convergence
plot(Blat.Iso.comb[[1]]$p, type = "l",
     ylim = c(-500,-250))
for (i in 2:100) {
  lines(Blat.Iso.comb[[i]]$p,
        col = rainbow(100)[i])
}

# load the first part of the order rates analysis (1-3)
load("../results/05a.order.analysis-1-3.RData")

# rename results accordingly
OrderRates1to3 <- results.orders

# rename Blattodea as Blattodea sensu stricto
# this means Blattodea in its original decription (without Isoptera)
names(OrderRates1to3)[1] <- "Blattodea sensu stricto"

# remove all but the rates
rm(list = ls()[c(-1,-12)])

# look at the convergence
#Blattodea sensu stricto
plot(OrderRates1to3$`Blattodea sensu stricto`[[1]]$p,
     type = "l",
     ylim = c(-300,-140))
for (i in 2:100) {
  lines(OrderRates1to3$`Blattodea sensu stricto`[[i]]$p,
        col = rainbow(100)[i])
}

# there is a single tree (tree number 72) which has a dip in in the
# likelihood below -160 between 600 - 800 generations.

for(i in 1:100){
  if(sum(OrderRates1to3$`Blattodea sensu stricto`[[i]]$p[600:800] < -160) > 10){
    print(i)
    plot(OrderRates1to3$`Blattodea sensu stricto`[[i]]$p,
         type = "l",
         ylim = c(-300,-140),
         main = i)
  }
}

# Isoptera
plot(OrderRates1to3$Isoptera[[1]]$p,
     type = "l",
     ylim = c(-160,-100))
for (i in 2:100) {
  lines(OrderRates1to3$Isoptera[[i]]$p,
        col = rainbow(100)[i])
}

# Mantodea
plot(OrderRates1to3$Mantodea[[1]]$p,
     type = "l",
     ylim = c(-110,-60))
for (i in 2:100) {
  lines(OrderRates1to3$Mantodea[[i]]$p,
        col = rainbow(100)[i])
}

# load the second part of the order rates analysis (4-5)
load("../results/05b.order.analysis-4-5.RData")

# rename results accordingly
OrderRates4to5 <- results.orders

# remove all but the rates
rm(list = ls()[c(-1,-12,-13)])

# look at the convergence
# orthoptera
plot(OrderRates4to5$Orthoptera[[1]]$p,
     type = "l",
     ylim = c(-90,-40))
for (i in 2:100) {
  lines(OrderRates4to5$Orthoptera[[i]]$p,
        col = rainbow(100)[i])
}

# phasmatodea
plot(OrderRates4to5$Phasmatodea[[1]]$p,
     type = "l",
     ylim = c(-180,-60))
for (i in 2:100) {
  lines(OrderRates4to5$Phasmatodea[[i]]$p,
        col = rainbow(100)[i])
}

# in phasmatodea there are three trees which shows odd behavior
for(i in 1:100){
  if(sum(OrderRates4to5$Phasmatodea[[i]]$p > -100) > 500){
    print(i)
    plot(OrderRates4to5$Phasmatodea[[i]]$p,
         type = "l",
         ylim = c(-120,-60),
         main = i)
  }
}

# there is allso two more tree which has not converged 
for(i in 1:100){
  if(var(OrderRates4to5$Phasmatodea[[i]]$p[-c(1:250)]) > 2){
    print(i)
    plot(OrderRates4to5$Phasmatodea[[i]]$p,
         type = "l",
         ylim = c(-180,-60),
         main = i)
  }
}

# these trees are 27,28,39 66 and 81.
# need rerun the MCMC for these trees before combining all the
# data into single file.
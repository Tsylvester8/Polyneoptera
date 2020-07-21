# Terrence Sylvester
# 1st April 2019

# Burnin

## read in data
# load the Blat-iso-combined rates
load("../results/05c.order.analysis-blat.iso.combined.RData")

# rename results accordingly
Blat.Iso.comb <- results

# remove all but the rates
rm(list = ls()[-1])

# load the first part of the order rates analysis (1-3)
load("../results/05a.order.analysis-1-3.RData")

# rename results accordingly
OrderRates1to3 <- results.orders

# rename Blattodea as Blattodea sensu stricto
# this means Blattodea in its original decription (without Isoptera)
names(OrderRates1to3)[1] <- "Blattodea sensu stricto"

# remove all but the rates
rm(list = ls()[c(-1,-12)])

# load the second part of the order rates analysis (4-5)
load("../results/07.order.rates.4-5.corrected.RData")

# rename results accordingly
OrderRates4to5 <- results.orders

# remove all but the rates
rm(list = ls()[c(-1,-2,-3)])

# make a single file with all the order rates which can be used for burnin
order.rates <- vector(mode = "list", length = 6)
names(order.rates) <- c("Blattodea", names(OrderRates1to3))

for (i in 1:6) {
  if(names(order.rates)[i] == "Blattodea"){
    order.rates[[i]] <- Blat.Iso.comb
  }
  else{
    hit <- which(names(order.rates)[i]==names(OrderRates1to3))
    order.rates[[i]] <- OrderRates1to3[[hit]]
  }
}

# remove all else but keep order rates
rm(list = ls()[-8])

# burnin 
burnin <- .25
iter <- 1000
burn <- -1:-(iter*burnin)


post.burnin <- vector("list", length=6)

names(post.burnin) <- c(names(order.rates))

for (i in 1:6) {
  x <- order.rates[[i]][[1]][burn,]
  for (j in 2:100) {
    x <- rbind(x, order.rates[[i]][[j]][burn,])
  }
  post.burnin[[i]] <- x
}

# remove all else but keep post.burnin
rm(list = ls()[-7])

# save image
save.image("../results/08.post.burnin.RData")

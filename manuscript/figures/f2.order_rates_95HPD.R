# Terrence Sylvester
# 4 February 2019
# Comparison of chromosome rates for the order level analysis

# load required libraries
library(coda)
library(viridis)

par(mfrow = c(1,3), mar = c(4.5,4.5,1,1))

# load the post burnin from the MCMC run
load(file = "../../analysis/results/08.post.burnin.RData")

# remove blatodea sensu stricto and isoptera
post.burnin <- post.burnin[c(-2,-3)]

# get colors
cols <- viridis(n = length(post.burnin))
cols2 <- viridis(n = length(post.burnin), alpha = .5)
#### comparison of fission rates across orders ####

# make the plot window
ymax <- 50
ymin <- -ymax * .1

plot(x = NA,
     y = NA,
     xlim = c(0,.85),
     ylim = c(ymin,ymax),
     xlab = "Rate of chromosome gains (per MYA)",
     ylab = "Density",
     main = "",
     cex.lab = 1)

counter <- ymin/length(post.burnin)

for(i in 1:length(post.burnin)){
  asc.density <- density(post.burnin[[i]]$asc1)
  if(i == 3){
    asc.density <- density(post.burnin[[i]]$asc1, bw = .0075)
  }
  polygon(asc.density,
          col = cols2[i],
          border = cols[i])
  HPD <- HPDinterval(as.mcmc(post.burnin[[i]]$asc1))
  segments(x0 = HPD[1], x1 = HPD[2],
           y0 = counter, y1 = counter,
           lwd = 4, col = cols[i])
  counter <- counter + (ymin/length(post.burnin))

}

#### comparison of fusion rates across orders ####

# make the plot window
ymax <- 30
ymin <- -ymax * .1

plot(x = NA, y = NA, xlim = c(0,1), ylim = c(ymin,ymax),
     xlab = "Rate of chromosome losses (per MYA)",
     ylab = "Density",
     main = "", cex.lab = 1)

counter <- ymin/length(post.burnin)

for(i in 1:length(post.burnin)){
  desc.density <- density(post.burnin[[i]]$desc1)
  if(i == 3){
    desc.density <- density(post.burnin[[i]]$desc1, bw = .011)
  }
  polygon(x = desc.density$x, y = desc.density$y,
          col = cols2[i], border = cols[i])
  HPD <- HPDinterval(as.mcmc(post.burnin[[i]]$desc1))
  segments(x0 = HPD[1], x1 = HPD[2],
           y0 = counter, y1 = counter,
           lwd = 4, col = cols[i])
  counter <- counter + (ymin/length(post.burnin))
}

#### comparision of polyploidy rates across orders ####

# make the plot window
ymax <- 30
ymin <- -ymax * .1

plot(x = NA,y = NA,xlim = c(0,1), ylim = c(ymin,ymax),
     xlab = "Rate of polyploidy (per MYA)",
     ylab = "Density",
     main = "", cex.lab = 1)
counter <- ymin/length(post.burnin)

for(i in 1:length(post.burnin)){
  poly.density <- density(post.burnin[[i]]$pol1)
  if(i == 1){
    poly.density <- density(post.burnin[[i]]$pol1, bw = .0135)
  }
  polygon(x = poly.density$x, y = poly.density$y,
          col = cols2[i], border = cols[i])
  HPD <- HPDinterval(as.mcmc(post.burnin[[i]]$pol1))
  segments(x0 = HPD[1], x1 = HPD[2],
           y0 = counter, y1 = counter,
           lwd = 4, col = cols[i])
  counter <- counter + (ymin/length(post.burnin))
}

### legend ###
y.counter <- 25
for (i in 1:length(post.burnin)) {
    points(x = .6, y = y.counter,
         pch = 15, col = cols[i],
         cex = 2)
    text(x = .6, y = y.counter,
       labels = names(post.burnin)[i], pos = 4,
       cex = 1)

y.counter <- y.counter - 2.5
}


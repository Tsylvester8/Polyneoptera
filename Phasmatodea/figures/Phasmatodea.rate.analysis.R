# Terrence Sylvester
# March 18 2019
# pradakshanas@gmail.com

# phasmatodea rates visualisation

# read in data
phas.rates <- read.csv(file = "../results/pruned.phas.rates.csv", as.is = T)

# load libraries
library(coda)

# state 1 is sexual reproduction
# state 2 is asexual reproduction

# differences of rates in each state (state 1 - state 2)
# chromosome gains / fissions
gain_difference <- phas.rates$asc1 - phas.rates$asc2
gain_HPD <- HPDinterval(as.mcmc(gain_difference))

# chromosome losses / fissions
loss_difference <- phas.rates$desc1 - phas.rates$desc2
loss_HPD <- HPDinterval(as.mcmc(loss_difference))

# polyploidy
poly_difference <- phas.rates$pol1 - phas.rates$pol2
poly_HPD <- HPDinterval(as.mcmc(poly_difference))

# devide the plot window
par(mfcol = c(1,3))

# plot of chromosome gains
plot(density(gain_difference),
     col = "blue",
     ylim = c(-1.5,15),
     xlim = c(-.15,.15),
     main = "Chromosome gains",
     xlab = "Rate difference (MYA)",
     cex.lab = 1.5,
     cex.main = 2)

polygon(density(gain_difference),
        border = NA,
        col = rgb(0,0,1,.2))

segments(x0 = gain_HPD[1], y0 = -0.75, x1 = gain_HPD[2], y1 = -0.75, col = "blue", lwd = 4)

abline(v = 0, col = "red", lty = 2, lwd = 3)

# put text
text(x = 0.12, y = 13.5, labels = "Faster in\nsexual", pos = 2, cex = 2)

text(x = -0.12, y = 13.5, labels = "Faster in\nasexual", pos = 4, cex = 2)

# plot of chromosome losses
plot(density(loss_difference),
     col = "blue",
     ylim = c(-1.7,17),
     xlim = c(-.15,.15),
     main = "Chromosome losses",
     xlab = "Rate difference (MYA)",
     cex.lab = 1.5,
     cex.main = 2)

polygon(density(loss_difference),
        border = NA,
        col = rgb(0,0,1,.2))

segments(x0 = loss_HPD[1], y0 = -0.85, x1 = loss_HPD[2], y1 = -0.85, col = "blue", lwd = 4)

abline(v = 0, col = "red", lty = 2, lwd = 3)

# put text
text(x = 0.12, y = 15.3, labels = "Faster in\nsexual", pos = 2, cex = 2)

text(x = -0.12, y = 15.3, labels = "Faster in\nasexual", pos = 4, cex = 2)


# plot if polyploidy
plot(density(poly_difference),
     col = "blue",
     ylim = c(-3.2,32),
     xlim = c(-.1,.05),
     main = "Polyploidy",
     xlab = "Rate difference (MYA)",
     cex.lab = 1.5,
     cex.main = 2)

polygon(density(poly_difference),
        border = NA,
        col = rgb(0,0,1,.2))

segments(x0 = poly_HPD[1], y0 = -1.6, x1 = poly_HPD[2], y1 = -1.6, col = "blue", lwd = 4)

abline(v = 0, col = "red", lty = 2, lwd = 3)

# put text
text(x = 0.04, y = 28.8, labels = "Faster in\nsexual", pos = 2, cex = 2)

text(x = -0.08, y = 28.8, labels = "Faster in\nasexual", pos = 4, cex = 2)

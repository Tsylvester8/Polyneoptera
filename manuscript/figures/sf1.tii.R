# Terrence Sylvester
# 29 January 2019
# pradakshanas@gmail.com
# figure for taxonomic instability index

# read in data
dat <- read.csv("../../analysis/data/taxonomic.instability.index/taxonomic.instability.index.csv",
                as.is = T)

# make the plot
plot(x = dat$Num,
     y = sort(dat$Index),
     xlab = "Taxon",
     ylab = "Taxonomic instability index",
     pch = 16,
     cex = .5)

# show the cut off point

abline(h = 4870,
       col = "red",
       lty = 3,
       lwd = 3)

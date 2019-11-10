# figure looking at the full dataset to see any patterns of
# chromosome evolution

# load required pachages
library(beeswarm)
library(viridis)

# read in data
chroms <- read.csv(file = "../data/chrom.data/chroms.csv",
                   as.is = T)

# filter out missing data
chroms <- chroms[chroms$SCS != "",]
chroms <- chroms[!is.na(chroms$hap),]

# Rename sex chromosome systems accordingly

# "complex XY|homomorphic" should be "XY-multi"
chroms[chroms$SCS == "complex XY|homomorphic",]$SCS <- "XY-multi" 
chroms[chroms$SCS == "complex XY",]$SCS <- "XY-multi" 

# "complex XO" should be "XY-multi"
chroms[chroms$SCS == "complex XO",]$SCS <- "XY-multi"

# "XY|homomorphic" should be "XY"
chroms[chroms$SCS == "XY|homomorphic",]$SCS <- "XY" 
chroms[chroms$SCS == "XO|XY",]$SCS <- "XY"
chroms[chroms$SCS == "homomorphic",]$SCS <- "XY"

# Include Isoptera under Blattodea
chroms[chroms$order == "Isoptera",]$order <- "Blattodea" 

# Exclude orders which do not have variation in sex chromsome systems
chroms.new <- data.frame()

for(i in 1:length(unique(chroms$order))){
  chroms.temp <- chroms[chroms$order == unique(chroms$order)[i], ]
  if(length(unique(chroms.temp$SCS)) > 1){
    chroms.new <- rbind(chroms.new, chroms.temp)
  }
  rm(chroms.temp)
}

# calculate the means of the haploid numbers for different SCSs in 
# different orders
chroms.means <- matrix(data = NA,
                       nrow = length(unique(chroms.new$order)),
                       ncol = 3)

rownames(chroms.means) <- c(unique(chroms.new$order))
colnames(chroms.means) <- c("XO", "XY", "XY-multi")

for(i in 1:length(unique(chroms.new$order))){
  for(j in 1:3){
    chroms.means[i,j] <-  mean(chroms$hap[chroms$order == unique(chroms.new$order)[i] & 
                                            chroms$SCS == colnames(chroms.means)[j]])
  }
}

# calculate the standatd errors of these means
# chroms.std.error <- matrix(data = NA,
#                            nrow = length(unique(chroms.new$order)),
#                            ncol = 3)

# rownames(chroms.std.error) <- c(unique(chroms.new$order))
# colnames(chroms.std.error) <- c("XO.se", "XY.se", "XY-multi.se")

# for(i in 1:length(unique(chroms.new$order))){
#   for(j in 1:3){
#     if(length(chroms$hap[chroms$order == unique(chroms.new$order)[i] & 
#                          chroms$SCS == colnames(chroms.means)[j]]) > 1){
#       chroms.std.error[i,j] <- chroms.means[i,j] / sqrt(length(chroms$hap[chroms$order == unique(chroms.new$order)[i] & 
#                                                                             chroms$SCS == colnames(chroms.means)[j]]))
#     }
#   }
#   
# }

# order them alphabetically
chroms.means <- chroms.means[order(rownames(chroms.means)),]

# chroms.std.error <- chroms.std.error[order(rownames(chroms.std.error)),]

# do the beeswarm plot

beeswarm(chroms.new$hap~chroms.new$SCS + chroms.new$order,
         corral = "gutter",
         cex = .9,
         pch = 16,
         xlab = "",
         ylab = "Haploid number",
         axes = F,
         col = viridis(n=3, alpha = 0.4,
                       direction = 1),
         spacing = 0.2,
         ylim = c(-1,50))


# counters
mean.x0 <- 0.6
mean.x1 <- 1.4

# se.upper.x0 <- 1
# se.upper.x1 <- 1
# 
# se.lower.x0 <- 1
# se.lower.x1 <- 1

for (i in 1:6) {
  for(j in 1:3){
    # plot means
    segments(x0 = mean.x0,
             y0 = chroms.means[i,j],
             x1 = mean.x1,
             y1 = chroms.means[i,j],
             col = "black",
             lwd = 3,
             lend = 2)
    # plot error bars
    # upper
    # segments(x0 = se.upper.x0,
    #          y0 = chroms.means[i,j],
    #          x1 = se.upper.x1,
    #          y1 = chroms.means[i,j]+chroms.std.error[i,j],
    #          col = "grey50",
    #          lty = 5)
    # horizontal bar
    # segments(x0 = se.upper.x0+0.2,
    #          y0 = chroms.means[i,j]+chroms.std.error[i,j],
    #          x1 = se.upper.x1-0.2,
    #          y1 = chroms.means[i,j]+chroms.std.error[i,j],
    #          col = "grey50",
    #          lty = 5)
    # lower
    # segments(x0 = se.lower.x0,
    #          y0 = chroms.means[i,j],
    #          x1 = se.lower.x1,
    #          y1 = chroms.means[i,j]-chroms.std.error[i,j],
    #          col = "grey50",
    #          lty = 5)
    #horizontal bar
    # segments(x0 = se.lower.x0 +0.2,
    #          y0 = chroms.means[i,j]-chroms.std.error[i,j],
    #          x1 = se.lower.x1 -0.2,
    #          y1 = chroms.means[i,j]-chroms.std.error[i,j],
    #          col = "grey50",
    #          lty = 5)
    
    # counters
    mean.x0 <- mean.x0 + 1
    mean.x1 <- mean.x1 + 1
    
    # se.upper.x0 <- se.upper.x0 + 1
    # se.upper.x1 <- se.upper.x1 + 1
    # 
    # se.lower.x0 <- se.lower.x0 + 1
    # se.lower.x1 <- se.lower.x1 + 1
    
  }
}

abline(v = (seq(from = 3.5, to = 18, by = 3)),
       col = "grey50",
       lty = 3)

axis(side = 1,
     at = c(0, seq(from = 2,  to = 18, by = 3), 19),
     labels = c(NA ,row.names(chroms.means), NA),
     tick = T,
     las = 1)

axis(side = 2,
     at = seq(from = 0,  to = 50, by = 10),
     tick = T)

left.s <- 16
for(i in 1:3){
  rect(xleft = left.s,
       ybottom = c(46,43,40)[i] ,
       xright = left.s+.5,
       ytop = c(48,45,42)[i],
       col = viridis(n = 3, alpha = 0.5,
                     direction = 1)[i])
  
  text(x = left.s+.5,
       y = c(47,44,41)[i],
       labels = sort(unique(chroms.new$SCS))[i],
       pos = 4,
       cex = .9)
  
}


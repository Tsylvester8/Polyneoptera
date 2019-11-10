# Terrence Sylvester
# February 18 2019
# analysis of means of the preliminary data set

# load libraries
library(car)

chroms <- read.csv(file = "../data/chrom.data/chroms.csv",
                   as.is = T)

# convert the sex chromosome systems in the original dataset
# so that there is only three types of sex chromosome systems

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
chroms$order[chroms$order %in% c("Grylloblata","Mantophasmatodea")] <- "Notoptera"

# get the orders names
orders <- unique(chroms$order)

# make a table with 3 columns for anova
# column 1 order name
# column 2 sex chromosome system
# column 3 haploid chromosome number
anov.mat <- matrix(data = NA,
                   nrow = nrow(chroms),
                   ncol = 3)

# label the column names
colnames(anov.mat) <- c("order", "SCS", "hap")

anov.mat <- as.data.frame(anov.mat)

counter <- 0

# fill the matrix
for(i in 1:length(orders)){
  n.dat <- sum(chroms$order==orders[i])
  counter <- counter + n.dat
  if(i == 1){
    anov.mat$order[1 :counter] <- chroms$order[chroms$order==orders[i]]
    anov.mat$SCS[1 :counter] <- chroms$SCS[chroms$order==orders[i]]
    anov.mat$hap[1 :counter] <- chroms$hap[chroms$order==orders[i]]  
  }
  
  anov.mat$order[(counter - n.dat + 1) :(counter)] <- chroms$order[chroms$order==orders[i]]
  anov.mat$SCS[(counter - n.dat + 1) :(counter)] <- chroms$SCS[chroms$order==orders[i]]
  anov.mat$hap[(counter - n.dat + 1) :(counter)] <- chroms$hap[chroms$order==orders[i]]
}

# test for variance in chromosome number in each order
leveneTest(anov.mat$hap~anov.mat$order)

# post hoc test

# Because a Levene test is simply an ANOVA conducted on sample variance (residuals) 
# instead of sample means, you can calculate the residuals manually, 
# then run the ANOVA with a TukeyHSD test as a post-hoc.

# get the meadians for each group
# store medians
med <- c()

# get medians for each group
for(i in 1:length(orders)){
med[i] <- median(anov.mat$hap[anov.mat$order == orders[i]])
}

names(med) <- orders

# get the residuals
# residuals
residuals <- c()

# get the absolute value of residuals for each order
for(i in 1:length(orders)){
  residuals <-c(residuals, abs(anov.mat$hap[anov.mat$order == orders[i]] - med[i]))
}

# make a new column to include the residuals
anov.mat$resid <- residuals

# do the levine's test to see if there is a difference when using residuals
# (there should not be any)
levene.dat.aov <- aov(anov.mat$resid~anov.mat$order)
# get the results of the levine's test
summary(levene.dat.aov)

post.hoc <- TukeyHSD(levene.dat.aov)

# get get the order pairs which shows significant difference in the
# variance in the chromosome number
results <- post.hoc$`anov.mat$order`[post.hoc$`anov.mat$order`[,4] < 0.05,]

# get the variences of the orders
for(i in 1:length(orders)){
  print(paste(orders[i],
              round(var(x = anov.mat$hap[anov.mat$order==orders[i]]),
                    2), 
              sep = " :"))
}



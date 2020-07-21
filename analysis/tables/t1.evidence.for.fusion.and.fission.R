# Terrence Sylvester
# March 4, 2020
# pradakshanas@gmail.com

# This will generate the table wich looks at whether chromosome number
# provides evidence for fusions and fissons in transitions between
# sex chromosome systems

# read in chromosome data
dat <- read.csv("../data/chrom.data/evidence.fusions.csv", as.is=T)

foo <- aggregate(dat$autosomes, by=list(as.character(dat$genus), dat$SCS), mean)


results <- data.frame(as.character(unique(foo$Group.1)),rep(NA,345),rep(NA,345),rep(NA,345),stringsAsFactors = F)
colnames(results) <- c("genus", "XO","XY","multi")
types <- c("XO", "XY","complex XY")
for(j in 1:3){
  cur.dat <- foo[foo$Group.2==types[j],]
  for(i in 1:nrow(results)){
    if(results$genus[i] %in% cur.dat$Group.1){
      hit <- which(cur.dat$Group.1 == results$genus[i])
      results[i,j+1] <- cur.dat$x[hit]
    }
  }
}
bad.rows <- c()
for(i in 1:nrow(results)){
  if(sum(is.na(results[i,]))>1){
    bad.rows <- c(bad.rows, i)
  }
}
results <- results[-bad.rows,]

gen.counts <- as.data.frame(table(dat$genus), stringsAsFactors = F)
gen.counts <- gen.counts[gen.counts$Var1 %in% results$genus,]

ssize <- c()
for(i in 1:nrow(results)){
  hit <- which(gen.counts$Var1 == results$genus[i])
  ssize[i] <- gen.counts$Freq[hit]
}

results <- data.frame(results, ssize)


# write results and edit in excell
write.csv(results,file="evid.fus.csv")

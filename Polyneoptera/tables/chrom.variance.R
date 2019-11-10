dat <- read.csv(file = "../data/chrom.data/chroms.csv",
                as.is = T)

dat$order[dat$order == "Isoptera"] <- "Blattodea"

prelim.data <- matrix(nrow = length(unique(dat$order)),
                      ncol = 4)
colnames(prelim.data) <- c("Order", 
                           "Records*", 
                           "Chromosome range (min, max)",
                           "Variance")

for(i in 1:length(unique(dat$order))){
  
  prelim.data[i,1] <- unique(dat$order)[i]
  dat.new <- dat[dat$order == unique(dat$order)[i],]
  dat.new <- dat.new[!is.na(dat.new$hap),]
  
  prelim.data[i,2] <- paste(nrow(dat[dat$order == unique(dat$order)[i],]),
                            " (",
                            nrow(dat.new),
                            ")",
                            sep = "") 
  prelim.data[i,3] <- paste(range(dat.new$hap)[1], 
                            "-" ,
                            range(dat.new$hap)[2],
                            sep = " ")
  prelim.data[i,4] <- round(var(dat.new$hap), digits = 2)
  if(is.na(prelim.data[i,4])){
    prelim.data[i,4] <- "-"
  }
}

write.csv(x = prelim.data, file = "../tables/chrom.variance.csv",
          row.names = F)

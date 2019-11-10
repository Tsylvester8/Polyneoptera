# This will generate the table

dat <- read.csv("../data/chrom.data/chroms.csv", as.is=T)
dat <- dat[,c(1,3,6,10)]
dat.clean <- dat[dat$SCS %in% unique(dat$SCS)[c(1,3,4,9)],]
results <- data.frame(unique(dat.clean$genus))
gen.count <- as.data.frame(table(dat$genus))
pos.gen <- as.character(gen.count[gen.count$Freq>1, 1])
var.gen <- c()
for(i in 1:length(pos.gen)){
  scs.pres <- unique(dat.clean$SCS[dat.clean$genus == pos.gen[i]])
  if(length(scs.pres)>1) var.gen[i] <- T
  if(length(scs.pres)<2) var.gen[i] <- F
}

proc.gen <- dat.clean[dat.clean$genus %in% pos.gen[var.gen], ]                     

# convert "complex XO" and "complex XY" to multi
proc.gen$SCS[proc.gen$SCS == "complex XO"] <- "multi"
proc.gen$SCS[proc.gen$SCS == "complex XY"] <- "multi"

# make the table 
gen.table <- matrix(ncol = 7,
                    nrow = length(unique(proc.gen$genus)))

colnames(gen.table) <- c("order","genera","XO","XY","multi","fusion","fission")

for(i in 1:length(unique(proc.gen$order))){
  proc.gen.temp <- proc.gen[proc.gen$order == unique(proc.gen$order)[i],]
  gen.table[sum(!is.na(gen.table[,1]))+1:length(unique(proc.gen.temp$genus)),1] <- rep(unique(proc.gen$order)[i], length(unique(proc.gen.temp$genus)))
  gen.table[sum(!is.na(gen.table[,2]))+1:length(unique(proc.gen.temp$genus)),2] <- unique(proc.gen.temp$genus)
}

gen.table <- as.data.frame(gen.table,
                           stringsAsFactors = F)

pos.scs <- unique(proc.gen$SCS)

for(i in 1:nrow(gen.table)){
  proc.gen.temp <- proc.gen[proc.gen$genus %in% gen.table$genera[i],]
  # for XO
  gen.table$XO[i] <- round(mean(proc.gen.temp$hap[proc.gen.temp$SCS == "XO"]),
                            digits = 1)
  # for XY
  gen.table$XY[i] <- round(mean(proc.gen.temp$hap[proc.gen.temp$SCS == "XY"]),
                            digits = 1)
  # for multi
  gen.table$multi[i] <- round(mean(proc.gen.temp$hap[proc.gen.temp$SCS == "multi"]),
                              digits = 1)
}

# looking at the gen table we see that there is one unknown genera and 
# one Phasmatodea genera for which the data is awailable for XO species.
# remove these two from the final table

gen.table <- gen.table[-c(4,18),]

# save the gen.table
write.csv(x = gen.table,
          file = "../tables/gen.table.csv", 
          row.names = F)

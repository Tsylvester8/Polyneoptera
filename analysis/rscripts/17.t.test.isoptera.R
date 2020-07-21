# Terrence Sylvester
# pradakshanas@gmail.com
# 8 April 2020
# Isoptera T test

# make a data table
dat <- matrix(data = NA,
              nrow = 100,
              ncol = 2)

# name columnds
colnames(dat) <- c("fusion", "fission")

# read in results
for(i in 1:100){
  

results <- readLines(paste("../asr.chrom.num/Isoptera/",
                           i,
                           "/results/expectations.txt",
                           sep = ""))

hit <- which(results == "#+++++++++++++++++++++++++++++")

dat[i,1] <- as.numeric(gsub(pattern = ".*: ", replacement = "", x = results[hit[1]-1]))
dat[i,2] <- as.numeric(gsub(pattern = ".*: ", replacement = "", x = results[hit[2]-1]))

}

# t test
t.test(dat[,1],dat[,2])
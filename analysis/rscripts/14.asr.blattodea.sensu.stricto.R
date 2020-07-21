# Terrence Sylvester
# 2 April 2019

# Ancestral state recomstruction using all the trees for blattodea sensu stricto
#### process ####

# this function allows for matching and 
# random choice when multiple records
source("02.helper.functions.R")

# this function will help in making the controll file
source("12.make.control.R")

# load ape for read in trees in nexus format
library(ape)

# read in data to get chromosome numbers
dat <- GetData(trees = "../data/trees/posterior.trees.nex", 
               data = "../data/chrom.data/chroms.csv")


# read in trees
trees <- read.nexus("../data/trees/posterior.trees.nex")

# read in the results from post processed results
load(file = "../results/post.burnin.RData")


# make the input files
for(j in 1:100){
  dat.pruned <- dat[dat$order == "Blattodea",]
  tree.pruned <- keep.tip(phy = trees[[j]],
                          tip = dat.pruned$species)
  
  # write the tree file in newick format
  if(!dir.exists(paths = paste("../asr.chrom.num/", "Blattodea_sensu_stricto","/",j, sep = ""))){
    dir.create(path = paste("../asr.chrom.num/", "Blattodea_sensu_stricto","/",j, sep = ""), 
               recursive = T)
  }
  
  write.tree(phy = tree.pruned, 
             file = paste("../asr.chrom.num/", "Blattodea_sensu_stricto", "/", j, "/", "Blattodea_sensu_stricto",".tre", sep = ""))
  
  # make the chromosome counts file
  counts <- matrix(nrow = nrow(dat.pruned),
                   ncol = 2)
  counts[,1] <- paste(">", dat.pruned$species, sep = "")
  counts[,2] <- dat.pruned$haploid
  
  write.table(counts, 
              file = paste("../asr.chrom.num/", "Blattodea_sensu_stricto","/", j, "/", "Blattodea_sensu_stricto",".count.txt", sep = ""),
              quote = F,
              sep = "\n",
              row.names = F,
              col.names = F)
  
  # chromosome range
  # this should be equal to what used in ChromPlus
  rng <- c(range(dat.pruned$haploid, na.rm = T)[1] - 2,
           range(dat.pruned$haploid, na.rm = T)[2] + 2)
  
  # number of simulations
  nsim <- 0
  
  make.control(mainType = "Run_Fix_Param",
               outDir = paste("../asr.chrom.num/", "Blattodea_sensu_stricto","/",j,"/","results", sep = ""),
               dataFile = paste("../asr.chrom.num/", "Blattodea_sensu_stricto","/", j, "/", "Blattodea_sensu_stricto",".count.txt", sep = ""),
               treeFile = paste("../asr.chrom.num/", "Blattodea_sensu_stricto", "/", j, "/", "Blattodea_sensu_stricto",".tre", sep = ""),
               maxChrNum = rng[2],
               minChrNum = rng[1],
               branchMul = 1,
               simulationsNum = nsim,
               pars = list(gainConstR = mean(post.burnin$`Blattodea sensu stricto`$asc1[((j*750)-(749)):(j*750)]),
                           lossConstR = mean(post.burnin$`Blattodea sensu stricto`$desc1[((j*750)-(749)):(j*750)]),
                           duplConstR = mean(post.burnin$`Blattodea sensu stricto`$pol1[((j*750)-(749)):(j*750)])),
               control.path = paste("../asr.chrom.num/", "Blattodea_sensu_stricto","/", j, "/", "Blattodea_sensu_stricto",".control.txt", sep = ""))
  
  # run ChromEvol
  system2(command = "../asr.chrom.num/chromEvol.exe",
          args = paste("../asr.chrom.num/", "Blattodea_sensu_stricto","/", j, "/", "Blattodea_sensu_stricto",".control.txt", sep = ""),
          wait = T)
}


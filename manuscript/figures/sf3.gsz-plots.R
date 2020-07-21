# Plots

# load libraries
library(ggplot2)
library(viridis)
library(ggpubr)
library(chromePlus)


# load data
load("../results/polyneoptera-gsz-chromnum-regressions.RData")
load("../results/gsz.vs.tip.rate.RData")

source("../rscripts/helper.functions.R")

chroms <- GetData(trees = "../data/Trees/posterior.trees.nex",
                  data = "../data/chrom.data/chroms.csv")

# set colours
group.colors <- c("Blattodea"   = viridis(7)[1],
                  "Dermaptera"    = viridis(7)[2],
                  "Embiidina" = viridis(7)[3],
                  "Isoptera"  = viridis(7)[4],
                  "Mantodea"   = viridis(7)[5],
                  "Orthoptera"  = viridis(7)[6],
                  "Phasmatodea"    = viridis(7)[7])

# get the order names

for(i in 1:nrow(dat)){
  dat$order[i] <- sample(chroms$order[which(chroms$species %in%
                                              dat$species[i])], 1)
}

# this will give the species names that are present in phylogenetically
# corrected lm
for(i in 1:nrow(gsize)){
  if(sum(gsize$Species[i] %in% gsize.sorted$Species) == 1){
    gsize$phylo[i] <- 16 # these are species that are aslo present in phylogenetically
                         # corrected data set
  }else{
    gsize$phylo[i] <- 17 # these are species that are not present in phylogenetically
                         # corrected data set
  }
}

# make plots
# make the plot
# genome size vs chromosome number
p1 <-  ggplot(gsize, aes(y=Chrom.num, x=Mbp)) +
  geom_point(aes(colour = Order), alpha=.7, size=2, shape = gsize$phylo, position = "jitter") +
  scale_colour_viridis_d() +
  theme_bw() +
  #theme(text=element_text(size=15, hjust=0.5, vjust=0.5)) +
  xlab("Genome size (Mbp)") +
  ylab("Haploid number") +
  annotate("text", x = c(7500, 7500), y = c(50, 47), hjust = 0, size=2,
           colour = c("blue","red"), label = c(paste("p-value",
              round(gsz.lm.sum$coefficients[2,4], 3)), paste("p-value",
              round(phylo.gsz.lm.sum$coefficients[2,4], 3)))) +
  geom_abline(intercept = phylo.gsz.lm.sum$coefficients[1],
              slope = phylo.gsz.lm.sum$coefficients[2],
              colour = "red", size = .8, linetype = "dashed") +
  geom_abline(intercept = lm(gsize$Chrom.num~gsize$Mbp)$coefficients[1],
              slope = lm(gsize$Chrom.num~gsize$Mbp)$coefficients[2],
              colour = "blue", size = .8) +
  labs(tag = "A")

# genome size vs tip rate
p2 <-  ggplot(dat, aes(y=meanRate, x=gsz)) +
  geom_point(aes(colour = order), alpha=.7, size=2, position = "jitter") +
  scale_colour_viridis_d() +
  theme_bw() +
  #theme(text=element_text(size=15, hjust=0.5, vjust=0.5)) +
  xlab("Genome size (Mbp)") + ylab("Mean tip rate") +
  geom_abline(intercept = phylo.fit.sum$coefficients[1],
              slope = phylo.fit.sum$coefficients[2],
              colour = "red", size = .8, linetype = "dashed") +
  geom_abline(intercept = fit.sum$coefficients[1],
              slope = fit.sum$coefficients[2],
              colour = "blue", size = .8) +
  annotate("text", x=rep(6100, 2), y=c(.33,.30),hjust = 0,
           colour = c("red","blue"), size=2,
           label = c(paste("p-value", round(fit.sum$coefficients[2,4], 3)),
                     paste("p-value",
                           round(phylo.fit.sum$coefficients[2,4], 3)))) +
  labs(tag = "B") +
  scale_colour_manual(values=group.colors)

p3 <- ggplot(dat, aes(y=abs(meanRate), x=gsz)) +
  geom_point(aes(colour = order), alpha=.7, size=2, position = "jitter") +
  scale_colour_viridis_d() +
  theme_bw() +
  #theme(text=element_text(size=15, hjust=0.5, vjust=0.5)) +
  xlab("Genome size (Mbp)") + ylab("Absolute mean tip rate") +
  geom_abline(intercept = phylo.fit.abs.sum$coefficients[1],
              slope = phylo.fit.abs.sum$coefficients[2],
              colour = "red", size = .8, linetype = "dashed") +
  geom_abline(intercept = fit.abs.sum$coefficients[1],
              slope = fit.abs.sum$coefficients[2],
              colour = "blue", size = .8) +
  annotate("text", x=rep(6100, 2), y=c(.33,.31),hjust = 0,
           colour = c("red","blue"), size=2,
           label = c(paste("p-value",
                           round(phylo.fit.abs.sum$coefficients[2,4], 3)),
                     paste("p-value",
                           round(fit.abs.sum$coefficients[2,4], 3)))) +
  labs(tag = "C") +
  scale_colour_manual(values=group.colors)

# plot the resutls
ggarrange(p1,p2,p3,
          common.legend = T,
          legend.grob = get_legend(p1),
          legend = "right",
          nrow = 1,
          ncol = 3)

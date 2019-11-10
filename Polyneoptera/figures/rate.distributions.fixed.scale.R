# figure looking at chromosome fusions, fissions and 
# rates of polyploidy in the 5 orders

# load required libraries
library(ggplot2)
library(gridExtra)
library(viridis)

# load data
load("../results/post.burnin.RData")

# to store plots
plots <-  list()

# ggplot will produce plots which have lables wich scientific notation
# for very small and very large numbers. Following code will remove
# scientific notation

options(scipen = 999)

#make and store plots
for(i in 1:length(post.burnin)){
  
  plots[[i]] <- ggplot(data = post.burnin[[i]][sample(1:75000, size = 10000),], 
                       mapping = aes(x =asc1,
                                     y =desc1,
                                     color = pol1)) + 
    geom_point(size = .3,
               alpha = 0.7,
               shape = 16,
               show.legend = T) +
    scale_colour_viridis(breaks = seq(0,1, by = .2), limits = c(0,1)) +
    labs(tag = LETTERS[i],
         subtitle = names(post.burnin)[i],
         color = "Rate of Polyploidy \n(per MY)") + 
    xlab("Ascending") +
    ylab("Descending") +
    xlim(c(0,1))+
    ylim(c(0,1))+
    theme_grey()
  
  names(plots)[i] <- names(post.burnin)[i]
}

# plot the resutls 
grid.arrange(grobs = c(plots[c(1:6)]),
             ncol = 3)

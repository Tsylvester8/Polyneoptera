# genome size impacts

gsz <- read.csv("../data/genome.sizes/gsz.orders.csv", as.is=T)
rates <- read.csv("../results/order.rates.post.burnin.csv", as.is=T)
orders <- unique(rates$order)
order.sizes <- gsz$GBP[gsz$group %in% orders]
names(order.sizes) <- gsz$group[gsz$group %in% orders]


gain <- loss <- pol <- c()
for(i in 1:length(order.sizes)){
  gain[i] <- mean(rates$asc[rates$order==names(order.sizes)[i]])
  loss[i] <- mean(rates$desc[rates$order==names(order.sizes)[i]])
  pol[i] <- mean(rates$pol[rates$order==names(order.sizes)[i]])
}
names(gain) <- names(loss) <- names(pol) <- names(order.sizes)

summary(lm(gain~order.sizes))
summary(lm(loss~order.sizes))
summary(lm(pol~order.sizes))

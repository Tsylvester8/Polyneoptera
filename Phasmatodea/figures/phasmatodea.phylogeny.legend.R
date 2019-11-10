library(viridis)

cols <- viridis(1001)

plot(x = NULL,
     y = NULL,
     xlim = c(0,10),
     ylim = c(0,10),
     axes = F,
     xlab = "",
     ylab = "")

add.color.bar(leg = 5,
              cols = cols,
              prompt = F,
              title = "",
              lims = NULL,
              lwd = 20,
              subtitle = "",
              x = 0,
              y = 5)

library(ggplot2)
library(gridExtra)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

# Plot variates by time
P <- c()
for (i in 2:10) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  plot((7988:8087)/4, data, type = "l")
  p <- data.frame(x = (7988:8087)/4, y = data) %>% ggplot(aes(x = x, y = y)) + geom_line() + geom_smooth() + ggtitle(header) + xlab("Time")
  P <- c(P, list(p))
}
do.call(grid.arrange, c(P, nrow = 3, ncol = 3))

P <- c()
for (i in 11:ncol(inflation.variates.data)) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  p <- data.frame(x = (7988:8087)/4, y = data) %>% ggplot(aes(x = x, y = y)) + geom_line() + geom_smooth() + ggtitle(header) + xlab("Time")
  P <- c(P, list(p))
}
do.call(grid.arrange, c(P, nrow = 3, ncol = 3))

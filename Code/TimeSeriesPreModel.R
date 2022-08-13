library(urca)
library(tseries)
library(GGally)
library(lavaan)
library(lattice)
library(RColorBrewer)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

# Plot variates by time
layout(matrix(1:9, byrow=TRUE, ncol=3,nrow=3))
par(mar=c(3, 3, 1.2, 0.5))
for (i in 2:10) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  plot((7988:8087)/4, data, type = "l")
  mtext(paste(header))
}
for (i in 11:ncol(inflation.variates.data)) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  plot((7988:8087)/4, data, type = "l")
  mtext(paste(header))
}
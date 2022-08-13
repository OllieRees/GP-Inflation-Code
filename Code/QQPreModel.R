library(urca)
library(tseries)
library(GGally)
library(lavaan)
library(lattice)
library(RColorBrewer)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

# QQ
layout(matrix(1:9, byrow=TRUE, ncol=3,nrow=3))
par(mar=c(3, 3, 1.2, 0.5))
for (i in 2:10) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  qqnorm(inflation.variates.data[, i], main = "")
  qqline(inflation.variates.data[, i], col = "steelblue", lwd = 2)
  mtext(paste("QQ of ", header))
}
for (i in 11:ncol(inflation.variates.data)) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  qqnorm(inflation.variates.data[, i], main = "")
  qqline(inflation.variates.data[, i], col = "steelblue", lwd = 2)
  mtext(paste("QQ of ", header))
}
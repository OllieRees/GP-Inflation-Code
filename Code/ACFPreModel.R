library(forecast)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

# ACF
layout(matrix(1:9, byrow=TRUE, ncol=3,nrow=3))
par(mar=c(3, 3, 1.2, 0.5))
for (i in 2:10) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  d <- auto.arima(data)$arma[6]
  if (d > 0) data <- diff(data, difference = d)
  autocorr <- acf(data, type = "correlation", plot = F, lag.max = 50)
  plot(autocorr)
  mtext(paste("ACF of ", header))
}
for (i in 11:ncol(inflation.variates.data)) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  d <- auto.arima(data)$arma[6]
  if (d > 0) data <- diff(data, difference = d)
  autocorr <- acf(data, type = "correlation", plot = F, lag.max = 50)
  plot(autocorr)
  mtext(paste("ACF of ", header))
}
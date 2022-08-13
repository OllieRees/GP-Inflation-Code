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

# Correlation matrix plot
ggcorr(inflation.variates.data[, 2:19])

# Covariance matrix
inflation.variates.cov <- cov(inflation.variates.data[, 2:18])
levelplot(inflation.variates.cov, col.regions = rainbow, main = "Covariance Heatmap of Non-Indexed Independent Vars.", pretty = T, contour = T,
          scales = list(y = list(at = 1:17, labels = colnames(inflation.variates.data[, 2:18])), 
                        x = list(at = 1:17, labels = 1:17)),
          ylab = "Variables", xlab = "Variable index (bottom to top)"
)

# Check for sationarity
for (i in 2:ncol(inflation.variates.data)) {
  #print(colnames(inflation.variates.data)[i])
  #for (lag in c(1, 4, 12, 25)) print(tseries::adf.test(inflation.variates.data[, i], k = lag))
}

# Check for cointegration between all explanatory variables
# Johansen test
for (i in 2:(ncol(inflation.variates.data) - 1)) {
  for (j in (i + 1):ncol(inflation.variates.data)) {
    #print(colnames(inflation.variates.data)[i])
    #print(colnames(inflation.variates.data)[j])
    summary(ca.jo(data.frame(inflation.variates.data[, i], inflation.variates.data[, j]), type="trace", K=2, ecdet="none", spec="longrun"))
  }
}

# ACF
layout(matrix(1:9, byrow=TRUE, ncol=3,nrow=3))
par(mar=c(3, 3, 1.2, 0.5))
for (i in 2:10) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  autocorr <- acf(data, type = "correlation", plot = F, lag.max = 50)
  plot(autocorr)
  mtext(paste("ACF of ", header))
}
for (i in 11:ncol(inflation.variates.data)) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  autocorr <- acf(data, type = "correlation", plot = F, lag.max = 50)
  plot(autocorr)
  mtext(paste("ACF of ", header))
}

# PACF
layout(matrix(1:9, byrow=TRUE, ncol=3,nrow=3))
par(mar=c(3, 3, 1.2, 0.5))
for (i in 2:10) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  autocorr <- pacf(data, type = "correlation", plot = F, lag.max = 50)
  plot(autocorr)
  mtext(paste("PACF of ", header))
}
for (i in 11:ncol(inflation.variates.data)) {
  header <- colnames(inflation.variates.data)[i]
  header <- gsub(".", " ", header, fixed = T)
  data <- inflation.variates.data[, i]
  autocorr <- pacf(data, type = "correlation", plot = F, lag.max = 50)
  plot(autocorr)
  mtext(paste("PACF of ", header))
}

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
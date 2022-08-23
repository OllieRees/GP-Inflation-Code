library(MASS)
library(highfrequency)
library(ggplot2)
library(tidyr)
library(lattice)
library(RColorBrewer)
library(latticeExtra)

set.seed(30642531)

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)
# Remove time, CPIH, Bank Rate, IB rate, 10Y Gilt, Exchange rate by CPI and Lab Cost
x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)] 
x$CPIHLag[2:100] <- inflation.variates.data[1:99, 2]
x$CPIHLag[1] <- 2.9 # From ONS
x <- data.matrix(x)
time.train <- (7988:8087/4)[41:100]
x.train <- x[41:100, ]
y.train <- inflation.variates.data$CPIH[41:100]
time.test <- (7988:8087/4)[1:40]
x.test <- x[1:40, ]
y.test <- inflation.variates.data$CPIH[1:40]

l <- res$x$l
kVar <- res$x$kVar
oDoF <- res$x$oDoF
iVars <- res$x$iVars
wMu <- res$x$wsMean
wVar <- res$x$wsVar

squared.exponential.kernel.point <- function(x1, x2) {
  kVar * exp(- sum( (x1 - x2)^2 ) / (2 * l^2))
}

squared.exponential.kernel <- function(X, Y) {
  x.size <- nrow(X)
  y.size <- nrow(Y)
  K <- matrix(0, nrow = x.size, ncol = y.size)
  for(i in 1:x.size) {
    for(j in 1:y.size) {
      K[i, j] <- squared.exponential.kernel.point(X[i, ], Y[j, ])
    }
  }
  K
}

se.partial.diff.point <- function(x, y, i, j, d) {
  x1 <- x[i]
  x2 <- y[j]
  x1d <- x[i, d]
  x2d <- y[j, d]
  se.k <- exp(-(sum( (x1 - x2)^2 ))/(2 * l^2))
  g2 <- -1/l^2 * (x1d - x2d)
  kVar * se.k * g2
}

se.partial.diff <- function(x, y) {
  x.size <- nrow(x)
  y.size <- nrow(y)
  D <- ncol(x)
  IK <- c()
  alpha <- rnorm(y.size, wMu, wVar)
  for (i in 1:x.size) {
    samevar <- c()
    for (d in 1:D) {
      samedim <- c()
      for (j in 1:y.size) {
        samedim[j] <- se.partial.diff.point(x, y, i, j, d)
        if (is.finite(samedim[j]) == F) {
          samedim[j] <- 0
        }
      }
      samevar[d] <- sum(samedim * alpha) 
    }
    IK <- rbind(IK, samevar)
  }  
  colnames(IK) <- colnames(x)
  rownames(IK) <- 1:x.size
  IK
}

cov.mat.func <- function(x, y) {
  outputNoiseVar <- diag(rchisq(nrow(x), oDoF), nrow = nrow(x), ncol = nrow(y))
  K <- squared.exponential.kernel(x, y)
  d <- se.partial.diff(x, y)
  iVars <- matrix(iVars, nrow = ncol(x), ncol = ncol(x), byrow = T)
  inputNoiseVar <- sqrt((d %*% iVars %*% t(d))^2) # NXxD * DxD * DxNY
  K + outputNoiseVar + inputNoiseVar
}

# Covariance heatmap
cov.prior <- cov.mat.func(x.test, x.test)

levelplot(cov.prior, col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Prior Covariance Heatmap", xlab = "", ylab = "")
levelplot(squared.exponential.kernel(x.test, x.test), col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Prior Kernel Heatmap (data)", xlab = "", ylab = "")
levelplot(squared.exponential.kernel(matrix(0:100), matrix(0:100)), col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Prior Kernel Heatmap (test)", xlab = "", ylab = "")

# plot prior
test.values <- data.frame(x = time.test, trueCPIH = y.test)
y.prior <- matrix(y.test, nrow = 40)
for (i in 1:10) {
  y.sample.prior <- MASS::mvrnorm(1, rep(0, 40), makePsd(cov.prior))
  test.values[gsub(" ", "", paste("CPIH", i))] <- y.sample.prior
}

test.values <- gather(test.values, variable, y, trueCPIH, CPIH1, CPIH2, CPIH3, CPIH4, CPIH5, CPIH6, CPIH7, CPIH8, CPIH8, CPIH9, CPIH10) 
ggplot(test.values, aes(x = x, y = y, color = variable)) + geom_line() + geom_smooth() + ggtitle("GP Prior") + xlab("Time") + ylab("CPIH") 

# Posterior
print(dim(squared.exponential.kernel(x.train, x.test)))
print(dim(solve(cov.prior)))
mean <- squared.exponential.kernel(x.train, x.test) %*% solve(cov.prior) %*% y.test
cov.post <- {
  a <- squared.exponential.kernel(x.train, x.train) 
  b <- squared.exponential.kernel(x.train, x.test)
  c <- solve(cov.prior)
  d <- squared.exponential.kernel(x.test, x.train)
  a - b %*% c %*% d
}

train.values <- data.frame(x = time.train, trueCPIH = y.train)
y.prior <- matrix(y.test, nrow = 40)
for (i in 1:10) {
  y.post <- MASS::mvrnorm(1, mean, makePsd(cov.post))
  train.values[gsub(" ", "", paste("CPIH", i))] <- y.post
}

train.values <- gather(train.values, variable, y, trueCPIH, CPIH1, CPIH2, CPIH3, CPIH4, CPIH5, CPIH6, CPIH7, CPIH8, CPIH8, CPIH9, CPIH10) 
ggplot(train.values, aes(x = x, y = y, color = variable)) + geom_line() + geom_smooth() + ggtitle("GP Posterior") + xlab("Time") + ylab("CPIH")
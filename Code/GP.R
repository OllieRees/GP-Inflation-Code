library(MASS)
library(highfrequency)
library(tidyr)
library(ggplot2)
library(lattice)
library(RColorBrewer)
library(latticeExtra)

set.seed(30642531)

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)
x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)] 
x$CPIHLag[2:100] <- inflation.variates.data[1:99, 2]
x$CPIHLag[1] <- 2.9
x <- data.matrix(x)
time.train <- (7988:8087/4)[41:100]
x.train <- x[41:100, ]
y.train <- inflation.variates.data$CPIH[41:100]
time.test <- (7988:8087/4)[1:40]
x.test <- x[1:40, ]
y.test <- inflation.variates.data$CPIH[1:40]

data.frame(x = time.test, y = y.test) %>% ggplot(aes(x = x, y = y)) + geom_line() + geom_smooth()

p <- res$x$p
lp <- res$x$lp
leq <- res$x$leq
kVar <- res$x$kVar
oDoF <- res$x$oDoF
iVars <- res$x$iVars
wMu <- res$x$wsMean
wVar <- res$x$wsVar

local.periodic.kernel.point <- function(x1, x2) {
  euclidean <- sum( (x1 - x2)^2 )
  periodic.sin.expr <- sin(pi/p * sqrt(euclidean))
  periodic.kernel <- exp(- 2/lp^2 * periodic.sin.expr^2)
  se.kernel <- exp(-euclidean / (2 * leq^2))
  kVar * periodic.kernel * se.kernel
}

local.periodic.kernel <- function(X, Y) {
  x.size <- nrow(X)
  y.size <- nrow(Y)
  K <- matrix(0, nrow = x.size, ncol = y.size)
  for(i in 1:x.size) {
    for(j in 1:y.size) {
      K[i, j] <- local.periodic.kernel.point(X[i, ], Y[j, ])
    }
  }
  K
}

local.periodic.partial.diff.point <- function(x, y, i, j, d) {
  x1 <- x[i]
  x2 <- y[j]
  x1d <- x[i, d]
  x2d <- y[j, d]
  p.sin <- sin(pi/p * sqrt(sum((x1 - x2)^2)) )
  p.k <- exp(-2/lp^2 * p.sin^2)
  se.k <- exp(-(sum( (x1 - x2)^2 ))/(2 * leq^2))
  g1 <- -4/lp^2 * p.sin * cos(pi/p * sqrt(sum((x1 - x2)^2))) * pi/p * (x1d - x2d)/sqrt(sum((x1 - x2)^2))
  g2 <- -1/leq^2 * (x1d - x2d)
  kVar * p.k * se.k * (g1 + g2)
}

local.periodic.partial.diff <- function(x, y) {
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
        samedim[j] <- local.periodic.partial.diff.point(x, y, i, j, d)
        if (is.finite(samedim[j]) == F) {
          samedim[j] <- 0
        }
      }
      samevar[d] <- sum(samedim * alpha)
    }
    IK <- rbind(IK, samevar)
  }  
  colnames(IK) <- colnames(x)
  rownames(IK) <- 1:(x.size)
  IK
}

cov.mat.func <- function(x, y) {
  outputNoiseVar <- diag(rchisq(nrow(x), oDoF), nrow = nrow(x), ncol = nrow(y))
  K <- local.periodic.kernel(x, y)
  d <- local.periodic.partial.diff(x, y)
  iVars <- matrix(iVars, nrow = ncol(x), ncol = ncol(x), byrow = T)
  print(d)
  inputNoiseVar <- d %*% iVars %*% t(d) # NXxD * DxD * DxNX
  K + outputNoiseVar + inputNoiseVar
}

cov.prior <- cov.mat.func(x.test, x.test)

# Covariance heatmap 
levelplot(cov.prior, col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Prior Covariance Heatmap", xlab = "", ylab = "")
levelplot(local.periodic.kernel(x.test, x.test), col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Prior Kernel Heatmap (data)", xlab = "", ylab = "")
levelplot(local.periodic.kernel(matrix(0:100), matrix(0:100)), col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Prior Kernel Heatmap (test)", xlab = "", ylab = "")

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
mean <- local.periodic.kernel(x.train, x.test) %*% solve(cov.prior) %*% y.test
cov.post <- {
  a <- local.periodic.kernel(x.train, x.train) 
  b <- local.periodic.kernel(x.train, x.test)
  c <- solve(cov.prior)
  d <- local.periodic.kernel(x.test, x.train)
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

# Mean and Variance of MAE, MSE, and RMSE
MAE.post <- c()
MSE.post <- c()
RMSE.post <- c()
for (i in 1:1000) {
  y.post <- MASS::mvrnorm(1, mean, makePsd(cov.post))
  r <- y.train - y.post
  MAE.post[i] <- mean(abs(r))
  MSE.post[i] <- mean((r)^2)
  RMSE.post[i] <- mean(sqrt((r)^2))
}

metrics.df <- data.frame(mean = c(mean(MAE.post), mean(MSE.post), mean(RMSE.post)), var = c(var(MAE.post), var(MSE.post), var(RMSE.post)))
rownames(metrics.df) <- c("MAE", "MSE", "RMSE")

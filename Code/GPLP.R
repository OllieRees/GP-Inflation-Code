library(MASS)
library(highfrequency)
library(tidyr)
library(ggplot2)
library(lattice)
library(RColorBrewer)
library(latticeExtra)
library(fracdiff)
library(forecast)

set.seed(30642531)

# Use optimal hyper-parameters flag
useOptimal = T

# Posterior data flag
trainingAsPostData = T

# Uses GP Allocation (last 60% as training data)
useGPAllocation = T

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)] 

x$CPIHLag[2:100] <- inflation.variates.data[1:99, 2]
x$CPIHLag[1] <- 2.9

x <- data.matrix(x)

inds.train <- 41:100
inds.test <- 1:40

if (useGPAllocation == F) {
  inds.train <- 1:80
  inds.test <- 81:100
}

if (trainingAsPostData) {
  time.train <- (7988:8087/4)[inds.train]
  x.train <- x[inds.train, ]
  y.train <- inflation.variates.data$CPIH[inds.train]

  time.test <- (7988:8087/4)[inds.test]
  x.test <- x[inds.test, ]
  y.test <- inflation.variates.data$CPIH[inds.test]
} else {
  time.train <- (7988:8087/4)[inds.test]
  x.train <- x[inds.test, ]
  y.train <- inflation.variates.data$CPIH[inds.test]
  
  time.test <- (7988:8087/4)[inds.train]
  x.test <- x[inds.train, ]
  y.test <- inflation.variates.data$CPIH[inds.train]
}

p <- res$x$p
lp <- res$x$lp * ifelse(useOptimal, 2, 1)
leq <- res$x$leq * ifelse(useOptimal, 2, 1)
kVar <- res$x$kVar 
oDoF <- res$x$oDoF # (or 0)

local.periodic.kernel.point <- function(x1, x2) {
  euclidean <- sum( (x1 - x2)^2 )
  periodic.sin.expr <- sin(pi/p * abs(x1 - x2))
  periodic.kernel <- exp(- 2/lp^2 * sum(periodic.sin.expr^2))
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

cov.mat.func <- function(x, y) {
  outputNoiseVar <- diag(rchisq(nrow(x), oDoF), nrow = nrow(x), ncol = nrow(y))
  K <- local.periodic.kernel(x, y)
  K + outputNoiseVar
}

cov.prior <- cov.mat.func(x.test, x.test)

if (useGPAllocation) {
  mean.prior <- rep(0, 60)
  if (trainingAsPostData) {
    mean.prior <- rep(0, 40)
  }
} else {
  mean.prior <- rep(0, 80)
  if (trainingAsPostData) {
    mean.prior <- rep(0, 20)
  }
}

# Covariance heatmap 
levelplot(cov.prior, 
          col.regions = rainbow, 
          scales = list(
            x = list(at = 0:10 * 4 + 1, labels = 1997:2007),
            y = list(at = 0:10 * 4 + 1, labels = 1997:2007)), 
          main = "Prior Covariance Heatmap (Locally Periodic)", xlab = "", ylab = "")
levelplot(local.periodic.kernel(x.test, x.test), 
          col.regions = rainbow, 
          scales = list(
            x = list(at = 0:10 * 4 + 1, labels = 1997:2007),
            y = list(at = 0:10 * 4 + 1, labels = 1997:2007)), 
          main = "Prior Covariance caused by Locally Periodic Kernel", xlab = "", ylab = "")
levelplot(local.periodic.kernel(matrix(0:100), matrix(0:100)), col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Locally Periodic Kernel ", xlab = "", ylab = "")

# plot Kernel
test.values <- data.frame(x = 1:100)
for (i in 1:10) {
  y.kern <- MASS::mvrnorm(1, rep(0, 100), local.periodic.kernel(matrix(1:100), matrix(1:100)))
  test.values[gsub(" ", "", paste("Draw", i))] <- y.kern
}

test.values <- gather(test.values, variable, y, Draw1, Draw2, Draw3, Draw4, Draw5, Draw6, Draw7, Draw8, Draw9, Draw10) 
ggplot(test.values, aes(x = x, y = y, color = variable)) + geom_line() + ggtitle("Locally Periodic Kernel Draw")

# plot prior
test.values <- data.frame(x = time.test, trueCPIH = y.test)
for (i in 1:10) {
  y.sample.prior <- MASS::mvrnorm(1, mean.prior, makePsd(cov.prior))
  test.values[gsub(" ", "", paste("CPIH", i))] <- y.sample.prior
}

test.values <- gather(test.values, variable, y, trueCPIH, CPIH1, CPIH2, CPIH3, CPIH4, CPIH5, CPIH6, CPIH7, CPIH8, CPIH8, CPIH9, CPIH10) 
ggplot(test.values, aes(x = x, y = y, color = variable)) + geom_line() + geom_smooth() + ggtitle("GP Prior (Locally Periodic)") + xlab("Time") + ylab("CPIH") 

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
for (i in 1:10) {
  y.post <- MASS::mvrnorm(1, mean, makePsd(cov.post))
  train.values[gsub(" ", "", paste("CPIH", i))] <- y.post
}

train.values <- gather(train.values, variable, y, trueCPIH, CPIH1, CPIH2, CPIH3, CPIH4, CPIH5, CPIH6, CPIH7, CPIH8, CPIH8, CPIH9, CPIH10) 
ggplot(train.values, aes(x = x, y = y, color = variable)) + geom_line() + geom_smooth() + ggtitle("GP Posterior (Locally Periodic)") + xlab("Time") + ylab("CPIH")

# Mean and Variance of MAE, MSE, and RMSE
MSE.post <- c()
MAPE.post <- c()
RMSE.post <- c()
for (i in 1:10000) {
  y.post <- MASS::mvrnorm(1, mean, makePsd(cov.post))
  r <- y.train - y.post
  MSE.post[i] <- mean((r)^2)
  MAPE.post[i] <- mean(abs(r/y.train))
  RMSE.post[i] <- mean(sqrt((r)^2))
}

metrics.df <- data.frame(mean = c(mean(MSE.post), mean(MAPE.post), mean(RMSE.post)), var = c(var(MSE.post), var(MAPE.post), var(RMSE.post)))
rownames(metrics.df) <- c("MSE", "MAPE", "RMSE")
print(metrics.df)

marginal.likelihood.log = function(y, cov.mat) {
  N <- length(y)
  if (class(y)[1] == "matrix") {
    N <- nrow(y)
  }  
  
  cov.mat.det <- det(cov.mat)
  
  # If statement used to prevent certain values from causing issues with invertibility
  if (all(is.finite(cov.mat)) && is.finite(cov.mat.det) && log(abs(cov.mat.det), 10) <= 100) {
    l <- -0.5 * (t(y) %*% solve(cov.mat) %*% y + log(abs(det(cov.mat))) + N * log(2 * pi))[1]
    if (is.finite(l)) {
      return(l)
    }
  }
  -99999
}

print(marginal.likelihood.log(y.test, cov.prior))
print(marginal.likelihood.log(y.train, cov.mat.func(x.train, x.train)))
library(MASS)
library(highfrequency)
library(ggplot2)
library(tidyr)
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

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)]  # Remove time, CPIH, Bank Rate, IB rate, 10Y Gilt, Exchange rate by CPI and Lab Cost

x$CPIHLag[2:100] <- inflation.variates.data[1:99, 2]
x$CPIHLag[1] <- 2.9 # From ONS

x <- data.matrix(x)
ivars.cov <- cov(x, method = "kendall")/10000

if (trainingAsPostData) {
  time.train <- (7988:8087/4)[41:100]
  x.train <- x[41:100, ]
  y.train <- inflation.variates.data$CPIH[41:100]

  time.test <- (7988:8087/4)[1:40]
  x.test <- x[1:40, ]
  y.test <- inflation.variates.data$CPIH[1:40]
} else {
  time.test <- (7988:8087/4)[41:100]
  x.test <- x[41:100, ]
  y.test <- inflation.variates.data$CPIH[41:100]
  
  time.train <- (7988:8087/4)[1:40]
  x.train <- x[1:40, ]
  y.train <- inflation.variates.data$CPIH[1:40]
}

l <- res$x$l * ifelse(useOptimal, 2, 1)
kVar <- res$x$kVar
oDoF <- res$x$oDoF

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

cov.mat.func <- function(x, y) {
  outputNoiseVar <- diag(rchisq(nrow(x), oDoF), nrow = nrow(x), ncol = nrow(y))
  K <- squared.exponential.kernel(x, y)
  K + outputNoiseVar
}

cov.prior <- cov.mat.func(x.test, x.test)

mean.prior <- rep(0, 60)
if (trainingAsPostData) {
  mean.prior <- rep(0, 40)
}

# Covariance heatmap 
levelplot(cov.prior, 
          col.regions = rainbow, 
          scales = list(
            x = list(at = 0:10 * 4 + 1, labels = 1997:2007),
            y = list(at = 0:10 * 4 + 1, labels = 1997:2007)), 
          main = "Prior Covariance Heatmap (Squared Exponential)", xlab = "", ylab = "")
levelplot(squared.exponential.kernel(x.test, x.test), 
          col.regions = rainbow, 
          scales = list(
            x = list(at = 0:10 * 4 + 1, labels = 1997:2007),
            y = list(at = 0:10 * 4 + 1, labels = 1997:2007)), 
          main = "Prior Covariance caused by Squared Exponential Kernel", xlab = "", ylab = "")
levelplot(squared.exponential.kernel(matrix(0:100), matrix(0:100)), col.regions = rainbow, scales = list(x = list(at = c(), labels = c()), y = list(at = c(), labels = c())), main = "Squared Exponential Kernel ", xlab = "", ylab = "")

# plot Kernel
test.values <- data.frame(x = 1:100)
for (i in 1:10) {
  y.kern <- MASS::mvrnorm(1, rep(0, 100), squared.exponential.kernel(matrix(1:100), matrix(1:100)))
  test.values[gsub(" ", "", paste("Draw", i))] <- y.kern
}

test.values <- gather(test.values, variable, y, Draw1, Draw2, Draw3, Draw4, Draw5, Draw6, Draw7, Draw8, Draw9, Draw10) 
ggplot(test.values, aes(x = x, y = y, color = variable)) + geom_line() + ggtitle("Squared Exponential Kernel Draw")

# plot prior
test.values <- data.frame(x = time.test, trueCPIH = y.test)
y.prior <- matrix(y.test, nrow = 40)
for (i in 1:10) {
  y.sample.prior <- MASS::mvrnorm(1, mean.prior, makePsd(cov.prior))
  test.values[gsub(" ", "", paste("CPIH", i))] <- y.sample.prior
}

test.values <- gather(test.values, variable, y, trueCPIH, CPIH1, CPIH2, CPIH3, CPIH4, CPIH5, CPIH6, CPIH7, CPIH8, CPIH8, CPIH9, CPIH10) 
ggplot(test.values, aes(x = x, y = y, color = variable)) + geom_line() + geom_smooth() + ggtitle("GP Prior (Squared Exponential)") + xlab("Time") + ylab("CPIH") 

# Posterior
mean <- squared.exponential.kernel(x.train, x.test) %*% solve(cov.prior) %*% y.test
cov.post <- {
  a <- squared.exponential.kernel(x.train, x.train) 
  b <- squared.exponential.kernel(x.train, x.test)
  c <- solve(cov.prior)
  d <- squared.exponential.kernel(x.test, x.train)
  a - b %*% c %*% d
}

train.values <- data.frame(x = time.train, trueCPIH = y.train)
for (i in 1:10) {
  y.post <- MASS::mvrnorm(1, mean, makePsd(cov.post))
  train.values[gsub(" ", "", paste("CPIH", i))] <- y.post
}

train.values <- gather(train.values, variable, y, trueCPIH, CPIH1, CPIH2, CPIH3, CPIH4, CPIH5, CPIH6, CPIH7, CPIH8, CPIH8, CPIH9, CPIH10) 
ggplot(train.values, aes(x = x, y = y, color = variable)) + geom_line() + geom_smooth() + ggtitle("GP Posterior (Squared Exponential)") + xlab("Time") + ylab("CPIH")

# Mean and Variance of MAE, MSE, and RMSE
MSE.post <- c()
MAE.post <- c()
MAPE.post <- c()
RMSE.post <- c()
for (i in 1:10000) {
  y.post <- MASS::mvrnorm(1, mean, makePsd(cov.post))
  r <- y.train - y.post
  MSE.post[i] <- mean((r)^2)
  MAE.post[i] <- mean(abs(r))
  MAPE.post[i] <- mean(abs(r/y.train))
  RMSE.post[i] <- mean(sqrt((r)^2))
}

metrics.df <- data.frame(mean = c(mean(MSE.post), mean(MAE.post), mean(MAPE.post), mean(RMSE.post)), var = c(var(MSE.post), var(MAE.post), var(MAPE.post), var(RMSE.post)))
rownames(metrics.df) <- c("MSE", "MAE", "MAPE", "RMSE")
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
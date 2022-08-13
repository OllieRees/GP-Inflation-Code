library(MASS)
# Kernels
local.periodic.kernel.point <- function(x1, x2, p, lp, leq, var) {
  euclidean <- function(a, b) sqrt(sum( (a-b)^2 ))
  periodic.kernel <- exp(- 2 * sum(sin(pi/p * abs(x1 - x2))^2 / lp^2))
  se.kernel <- exp(euclidean(x1, x2) / (2 * leq^2))
  var * periodic.kernel * se.kernel
}

local.periodic.kernel <- function(X, Y, p, lp, leq, var) {
  x.size <- nrow(X)
  K <- matrix(0, nrow = x.size, ncol = x.size)
  for(i in 1:x.size) {
    for(j in 1:x.size) {
      K[i, j] <- local.periodic.kernel.point(X[i, ], Y[j, ], p, lp, leq, var)
    }
  }
  K
}

# Create Covariance matrix
cov.mat.func <- function(x, y, params) {
  outputNoiseVar <- diag(rep(params$oVar, nrow(x)))
  inputNoiseVar <- matrix(params$iVars, nrow = 100, ncol = 100, byrow = T) 
  K <- local.periodic.kernel(x, y, params$p, params$lp, params$leq, params$kVar)
  K + outputNoiseVar + inputNoiseVar
}

cov.prior <- cov.mat.func(x.test, y.test, params)

# Prior
y.prior <- MASS:rnormmvn(nrow(x.test), 0, cov.prior)

# Posterior
mean <- local.periodic.kernel(x.train, x.test, params$p, params$lp, params$leq, params$kVar) %*% solve(cov.prior) %*% y.prior
cov.post <- {
  a <- local.periodic.kernel(x.test, x.test, params$p, params$lp, params$leq, params$kVar) 
  b <- local.periodic.kernel(x.train, x.test, params$p, params$lp, params$leq, params$kVar)
  c <- solve(cov.prior)
  d <- local.periodic.kernel(x.test, x.train, params$p, params$lp, params$leq, params$kVar)
  a - b %*% c %*% d
}

y.post <- MASS:rnormmvn(nrow(x.train), mean, cov.post)
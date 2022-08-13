library(MASS)
library(ggplot2)
library(mlrMBO)
library(parallelMap)

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)
# Remove time, CPIH, Bank Rate, IB rate, 10Y Gilt, Exchange rate by CPI and Lab Cost
x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)] 
x$CPIHLag[2:100] <- inflation.variates.data[1:99, 2]
x$CPIHLag[1] <- 2.9 # From ONS
x <- data.matrix(x)
x <- x[41:100, ]
y <- inflation.variates.data$CPIH[41:100]

N <- nrow(x)
D <- ncol(x)

squared.exponential.kernel.point <- function(x1, x2, l, var) {
  var * exp(sum( (x1 - x2)^2 ) / (2 * l^2))
}

squared.exponential.kernel <- function(X, l, var) {
  x.size <- nrow(X)
  K <- matrix(0, nrow = x.size, ncol = x.size)
  for(i in 1:x.size) {
    for(j in 1:x.size) {
      K[i, j] <- squared.exponential.kernel.point(X[i, ], X[j, ], l, var)
    }
  }
  K
}

se.partial.diff.point <- function(x, i, j, d, l, var) {
  x1 <- x[i]
  x2 <- x[j]
  x1d <- x[i, d]
  x2d <- x[j, d]
  se.k <- exp(-(sum( (x1 - x2)^2 ))/(2 * l^2))
  g2 <- -1/l^2 * (x1d - x2d)
  var * se.k * g2
}

se.partial.diff <- function(x, l, var, alpha) {
  N <- nrow(x)
  D <- ncol(x)
  IK <- c()
  for (i in 1:N) {
    samevar <- c()
    for (d in 1:D) {
      samedim <- c()
      for (j in 1:N) {
        samedim[j] <- se.partial.diff.point(x, i, j, d, l, var)
        if (is.finite(samedim[j]) == F) {
          samedim[j] <- 0
        }
      }
      samevar[d] <- sum(samedim * alpha) 
    }
    IK <- rbind(IK, samevar)
  }  
  colnames(IK) <- colnames(x)
  rownames(IK) <- 1:nrow(x)
  IK
}


# Create Covariance matrix
cov.mat.func <- function(x, l, kVar, oVar, iVars, ws) {
  iVars <- matrix(iVars, nrow = ncol(x), ncol = ncol(x), byrow = T)
  outputNoiseVar <- diag(rep(oVar, nrow(x)))
  K <- squared.exponential.kernel(x, l, kVar)
  d <- se.partial.diff(x, l, kVar, ws)
  inputNoiseVar <- d %*% iVars %*% t(d) # NxD * DxD * DxN
  K + outputNoiseVar + inputNoiseVar
}

# Surrogate function: Log marginal likelihood of the score on the hyperparameters
marginal.likelihood.log = function(y, cov.mat) {
  N <- length(y)
  if (class(y)[1] == "matrix") {
    N <- nrow(y)
  }  

  cov.mat.det <- det(cov.mat)
  print(cov.mat.det)
  # If statement used to prevent certain values from causing issues with invertibility
  if (all(is.finite(cov.mat)) && is.finite(cov.mat.det) && log(abs(cov.mat.det), 10) <= 100) {
    l <- -0.5 * (t(y) %*% solve(cov.mat) %*% y + log(abs(det(cov.mat))) + N * log(2 * pi))[1]
    if (is.finite(l)) {
      return(l)
    }
  }
  -Inf
}

# Bayes optimisation using EI
obj.fun = makeSingleObjectiveFunction(
  name = "Marginal Log Likelihood",
  description = "Maximises the Marginal Log Likelihood on the output",
  
  minimize = F,
  #has.simple.signature = F,
  
  fn = function(ps) {
    l <- ps[1]
    kVar <- ps[2]
    oVar <- ps[3]
    iVars <- ps[4:(D^2 + 3)]
    ws <- ps[(D^2 + 4):length(ps)]
    cov.mat <- cov.mat.func(x, l, kVar, oVar, iVars, ws)
    marginal.likelihood.log(y, cov.mat)
  },
  
  par.set = makeParamSet(
    makeNumericParam(id = "l", lower = 0, upper = 100),
    makeNumericParam(id = "kVar", lower = 0, upper = 10),
    makeNumericParam(id = "oVar", lower = 0, upper = 10),
    makeNumericVectorParam(id = "iVars", ncol(x)^2, lower = 0, upper = 10),
    makeNumericVectorParam(id = "ws", nrow(x), lower = 0, upper = 1)
  )
)


#cov.mat <- cov.mat.func(x, rnorm(1, pi), rnorm(1, 10), rnorm(1, 10), rnorm(1), rnorm(1), rnorm(13^2), rnorm(20))
#z <- marginal.likelihood.log(y, cov.mat)

ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 100)
ctrl <- setMBOControlInfill(ctrl, opt.focussearch.points = 500, crit = makeMBOInfillCritEI()) 

parallelStartMulticore(cpus = 20, show.info = TRUE)

des <- generateDesign(n = (5 + D^2 + N) * 3, getParamSet(obj.fun), fun = lhs::geneticLHS)
des$y <- apply(des, 1, function(xs) obj.fun(list(xs[1], xs[2], xs[3], xs[4:(D^2 + 3)], xs[(D^2 + 4):ncol(des)])))

surr <- makeLearner("regr.km", predict.type = "se")

res <- mbo(obj.fun, design = des, control = ctrl, learner = surr)

parallelStop()

print(res)
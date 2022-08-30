library(MASS)
library(ggplot2)
library(mlrMBO)
library(parallelMap)
library(fracdiff)
library(forecast)

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

#for (column in c(4, 11)) {
#   data <- inflation.variates.data[, column]
#   d <- auto.arima(data)$arma[6]
#   inflation.variates.data[, column] <- diffseries(inflation.variates.data[, column], d)
#}

x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)]# Remove time, CPIH, Bank Rate, IB rate, 10Y Gilt, Exchange rate by CPI and Lab Cost
#x <- data.frame()
x$CPIHLag[2:100] <- inflation.variates.data[1:99, 2]
x$CPIHLag[1] <- 2.9 # From ONS

x <- data.matrix(x)
ivars.cov <- cov(x, method = "kendall")/10000

x <- x[41:100, ]
y <- inflation.variates.data$CPIH[41:100]


local.periodic.kernel.point <- function(x1, x2, p, lp, leq, var) {
  euclidean <- sum( (x1 - x2)^2 )
  periodic.sin.expr <- sin(pi/p * abs(x1 - x2))
  periodic.kernel <- exp(- 2/lp^2 * sum(periodic.sin.expr^2))
  se.kernel <- exp(-euclidean / (2 * leq^2))
  var * periodic.kernel * se.kernel
}

local.periodic.kernel <- function(X, p, lp, leq, var) {
  x.size <- nrow(X)
  K <- matrix(0, nrow = x.size, ncol = x.size)
  for(i in 1:x.size) {
    for(j in 1:x.size) {
      K[i, j] <- local.periodic.kernel.point(X[i, ], X[j, ], p, lp, leq, var)
    }
  }
  K
}

local.periodic.partial.diff.point <- function(x, i, j, d, p, lp, leq, var) {
  x1 <- x[i]
  x2 <- x[j]
  x1d <- x[i, d]
  x2d <- x[j, d]
  p.sin <- sin(pi/p * abs(x1 - x2))
  p.k <- exp(-2/lp^2 * sum(p.sin^2))
  se.k <- exp(-(sum( (x1 - x2)^2 ))/(2 * leq^2))
  g1 <- -4/lp^2 * sin(pi/p * abs(x1d - x2d)) * cos(pi/p * abs(x1d - x2d)) * pi/p * (x1d - x2d)/abs(x1d - x2d)
  g2 <- -1/leq^2 * (x1d - x2d)
  var * p.k * se.k * (g1 + g2)
}

local.periodic.partial.diff <- function(x, p, lp, leq, var, wMu, wVar) {
  N <- nrow(x)
  D <- ncol(x)
  IK <- c()
  alpha <- rnorm(N, wMu, wVar)
  for (i in 1:N) {
    samevar <- c()
    for (d in 1:D) {
      samedim <- c()
      for (j in 1:N) {
        samedim[j] <- local.periodic.partial.diff.point(x, i, j, d, p, lp, leq, var)
        if (is.finite(samedim[j]) == F) {
          samedim[j] <- 0
        }
      }
      samevar[d] <- sum(samedim * alpha)
    }
    IK <- rbind(IK, samevar)
  }  
  colnames(IK) <- colnames(x)
  rownames(IK) <- 1:N
  IK
}

# Create Covariance matrix
cov.mat.func <- function(x, p, lp, leq, kVar, oDof) {
  outputNoiseVar <- diag(rchisq(nrow(x), oDof))
  K <- local.periodic.kernel(x, p, lp, leq, kVar)
  #d <- local.periodic.partial.diff(x, p, lp, leq, kVar, wsMean, wsVar)
  #inputNoiseVar <- sqrt((d %*% ivars.cov %*% t(d))^2) # NxD * DxD * DxN
  K + outputNoiseVar
}

# Surrogate function: Log marginal likelihood of the score on the hyperparameters
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
  
# Bayes optimisation using EI
obj.fun = makeSingleObjectiveFunction(
    name = "Marginal Log Likelihood",
    description = "Maximises the Marginal Log Likelihood on the output",
    
    minimize = F,
    
    fn = function(ps) {
      p <- ps[1]
      lp <- ps[2]
      leq <- ps[3]
      kVar <- ps[4]
      oDoF <- ps[5]
      cov.mat <- cov.mat.func(x, p, lp, leq, kVar, oDoF)
      marginal.likelihood.log(y, cov.mat)
    },
    
    par.set = makeParamSet(
      makeNumericParam(id = "p", lower = 0, upper = 22),
      makeNumericParam(id = "lp", lower = 0, upper = 12),
      makeNumericParam(id = "leq", lower = 0, upper = 60),
      makeNumericParam(id = "kVar", lower = 0, upper = 3),
      
      makeNumericParam(id = "oDoF", lower = 0, upper = 2)
    )
)

parallelStartMulticore(cpus = 20, show.info = F)
  
ctrl <- makeMBOControl()
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())
ctrl <- setMBOControlTermination(ctrl, iters = 500)

des <- generateDesign(n = 25, getParamSet(obj.fun), fun = lhs::geneticLHS)
des$y <- apply(des, 1, function(xs) obj.fun(list(xs[1], xs[2], xs[3], xs[4], xs[5])))
  
configureMlr(on.learner.warning = "quiet", show.learner.output = FALSE)
surr <- makeLearner("regr.km", predict.type = "se")
  
res <- mbo(obj.fun, design = des, control = ctrl, learner = surr)

parallelStop()
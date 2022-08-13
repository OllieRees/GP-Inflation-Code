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

local.periodic.kernel.point <- function(x1, x2, p, lp, leq, var) {
  euclidean <- sum( (x1 - x2)^2 )
  periodic.sin.expr <- sin(pi/p * abs(x1 - x2))^2
  periodic.kernel <- exp(- 2/lp^2 * sum(periodic.sin.expr))
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

local.periodic.partial.diff <- function(x, p, lp, leq, var, alpha) {
  N <- nrow(x)
  D <- ncol(x)
  IK <- c()
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
cov.mat.func <- function(x, p, lp, leq, kVar, oVar, iVars, ws) {
  iVars <- matrix(iVars, nrow = ncol(x), ncol = ncol(x), byrow = T)
  outputNoiseVar <- diag(rep(oVar, nrow(x)))
  K <- local.periodic.kernel(x, p, lp, leq, kVar)
  d <- local.periodic.partial.diff(x, p, lp, leq, kVar, ws)
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

  # If statement used to prevent certain values from causing issues with invertibility
  if (all(is.finite(cov.mat)) && is.finite(cov.mat.det) && log(abs(cov.mat.det), 10) <= 100) {
    l <- -0.5 * (t(y) %*% solve(cov.mat) %*% y + log(abs(det(cov.mat))) + N * log(2 * pi))[1]
    if (is.finite(l)) {
      return(l)
    }
  }
  0
}

# Bayes optimisation using EI
obj.fun = makeSingleObjectiveFunction(
  name = "Marginal Log Likelihood",
  description = "Maximises the Marginal Log Likelihood on the output",
  
  minimize = F,
  #has.simple.signature = F,
  
  fn = function(ps) {
    p <- ps[1]
    lp <- ps[2]
    leq <- ps[3]
    kVar <- ps[4]
    oVar <- ps[5]
    iVars <- ps[6:(D^2 + 5)]
    ws <- ps[(D^2 + 6):length(ps)]
    cov.mat <- cov.mat.func(x, p, lp, leq, kVar, oVar, iVars, ws)
    marginal.likelihood.log(y, cov.mat)
  },
  
  par.set = makeParamSet(
    makeNumericParam(id = "p", lower = 0, upper = pi),
    makeNumericParam(id = "lp", lower = 0, upper = 100),
    makeNumericParam(id = "leq", lower = 0, upper = 100),
    makeNumericParam(id = "kVar", lower = 0, upper = 10),
    makeNumericParam(id = "oVar", lower = 0, upper = 10),
    makeNumericVectorParam(id = "iVars", ncol(x)^2, lower = 0, upper = 10),
    makeNumericVectorParam(id = "ws", nrow(x), lower = 0, upper = 1)
  )
)

ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 100)
ctrl <- setMBOControlInfill(ctrl, opt.focussearch.points = 500, crit = makeMBOInfillCritEI()) 

parallelStartMulticore(cpus = 20, show.info = TRUE)

des <- generateDesign(n = (5 + D^2 + N) * 3, getParamSet(obj.fun), fun = lhs::geneticLHS)
des$y <- apply(des, 1, function(xs) obj.fun(list(xs[1], xs[2], xs[3], xs[4], xs[5], xs[6:(D^2 + 5)], xs[(D^2 + 6):ncol(des)])))

surr <- makeLearner("regr.km", predict.type = "se")

res <- mbo(obj.fun, design = des, control = ctrl, learner = surr)

parallelStop()

print(res)
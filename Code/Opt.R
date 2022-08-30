library(MASS)
library(ggplot2)
library(mlrMBO)
library(parallelMap)
library(fracdiff)
library(forecast)

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)]# Remove time, CPIH, Bank Rate, IB rate, 10Y Gilt, Exchange rate by CPI and Lab Cost

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

# Create Covariance matrix
cov.mat.func <- function(x, p, lp, leq, kVar, oDof) {
  outputNoiseVar <- diag(rchisq(nrow(x), oDof))
  K <- local.periodic.kernel(x, p, lp, leq, kVar)
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
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

squared.exponential.kernel.point <- function(x1, x2, l, var) {
  var * exp(- sum( (x1 - x2)^2 ) / (2 * l^2))
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

# Create Covariance matrix
cov.mat.func <- function(x, l, kVar, oDoF) {
  outputNoiseVar <- diag(rchisq(nrow(x), oDoF))
  K <- squared.exponential.kernel(x, l, kVar)
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
  if (all(is.finite(cov.mat)) && is.finite(cov.mat.det) && log(abs(cov.mat.det), 10) <= 90) {
    l <- -0.5 * (t(y) %*% solve(cov.mat) %*% y + log(abs(det(cov.mat))) + N * log(2 * pi))[1]
    if (is.finite(l)) {
      return(l)
    }
  }
  -999999
}

# Bayes optimisation using EI
obj.fun = makeSingleObjectiveFunction(
  name = "Marginal Log Likelihood",
  description = "Maximises the Marginal Log Likelihood on the output",
  
  minimize = F,
  
  fn = function(ps) {
    l <- ps[1]
    kVar <- ps[2]
    oDoF <- ps[3]
    cov.mat <- cov.mat.func(x, l, kVar, oDoF)
    marginal.likelihood.log(y, cov.mat)
  },
  
  par.set = makeParamSet(
    makeNumericParam(id = "l", lower = 0, upper = 60),
    makeNumericParam(id = "kVar", lower = 0, upper = 3),
    makeNumericParam(id = "oDoF", lower = 0, upper = 2)
  )
)

parallelStartMulticore(cpus = 20, show.info = TRUE)

ctrl <- makeMBOControl()
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())
ctrl <- setMBOControlTermination(ctrl, iters = 500)

des <- generateDesign(n = 5 * 3, getParamSet(obj.fun), fun = lhs::geneticLHS)
des$y <- apply(des, 1, function(xs) obj.fun(list(xs[1], xs[2], xs[3])))

configureMlr(on.learner.warning = "quiet", show.learner.output = FALSE)
surr <- makeLearner("regr.km", predict.type = "se")

res <- mbo(obj.fun, design = des, control = ctrl, learner = surr)

parallelStop()
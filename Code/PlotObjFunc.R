library(ggplot2)

# Data
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)
# Remove time, CPIH, Bank Rate, IB rate, 10Y Gilt, Exchange rate by CPI and Lab Cost
x <- inflation.variates.data[, -c(1, 2, 12, 14, 15, 17, 18)] 
x$CPIHLag[2:100] <- inflation.variates.data[1:99, 2]
x$CPIHLag[1] <- 2.9 # From ONS
x <- data.matrix(x)
x <- x[1:50, ]
y <- inflation.variates.data$CPIH[1:50]
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
cov.mat.func <- function(x, p, lp, leq, kVar, oDof, iVars, wsMean, wsVar) {
  outputNoiseVar <- diag(rchisq(N, oDof))
  K <- local.periodic.kernel(x, p, lp, leq, kVar)
  d <- local.periodic.partial.diff(x, p, lp, leq, kVar, wsMean, wsVar)
  iVars <- matrix(iVars, nrow = ncol(x), ncol = ncol(x), byrow = T)
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
  -99999
}

opt.params <- res$x
objs <- c()

# p is the variable
ps <- (1:(pi * 100))/100
for (p in ps) {
  cov.mat <- cov.mat.func(x, p, opt.params$lp, opt.params$leq, opt.params$kVar, opt.params$oDoF, opt.params$iVars, opt.params$wsMean, opt.params$wsVar)
  objs <- append(objs, marginal.likelihood.log(y, cov.mat))
}
df <- data.frame(Period = ps, maringalLogLikelihood = objs)
ggplot(df, aes(x = Period, y = maringalLogLikelihood)) + geom_line(color = "red") + geom_smooth()

# lp is the variable
lps <- (1:(100 * 100))/100
for (lp in lps) {
  cov.mat <- cov.mat.func(x, opt.params$p, lp, opt.params$leq, opt.params$kVar, opt.params$oDoF, opt.params$iVars, opt.params$wsMean, opt.params$wsVar)
  objs <- append(objs, marginal.likelihood.log(y, cov.mat))
}
df <- data.frame(LengthScaleP = lps, MaringalLogLikelihood = objs)
ggplot(df, aes(x = LengthScaleP, y = MaringalLogLikelihood)) + geom_line(color = "red") + geom_smooth()

# leq is the variable
leqs <- (1:(100 * 100))/100
for (leq in leqs) {
  cov.mat <- cov.mat.func(x, opt.params$p, opt.params$lp, leq, opt.params$kVar, opt.params$oDoF, opt.params$iVars, opt.params$wsMean, opt.params$wsVar)
  objs <- append(objs, marginal.likelihood.log(y, cov.mat))
}
df <- data.frame(LengthScaleSE = leqs, maringalLogLikelihood = objs)
ggplot(df, aes(x = LengthScaleSE, y = maringalLogLikelihood)) + geom_line(color = "red") + geom_smooth()

# kVar is the variable
kVars <- (1:(10 * 30))/30
for (kVar in kVars) {
  cov.mat <- cov.mat.func(x, opt.params$p, opt.params$lp, opt.params$leq, kVar, opt.params$oDoF, opt.params$iVars, opt.params$wsMean, opt.params$wsVar)
  objs <- append(objs, marginal.likelihood.log(y, cov.mat))
}
df <- data.frame(KernelVariance = kVars, MaringalLogLikelihood = objs)
ggplot(df, aes(x = KernelVariance, y = MaringalLogLikelihood)) + geom_line(color = "red") + geom_smooth()

# oDoF is the variable
oDofs <- (1:(10 * 30))/30
for (oDoF in oDofs) {
  cov.mat <- cov.mat.func(x, opt.params$p, opt.params$lp, opt.params$leq, opt.params$kVar, oDoF, opt.params$iVars, opt.params$wsMean, opt.params$wsVar)
  objs <- append(objs, marginal.likelihood.log(y, cov.mat))
}
df <- data.frame(OutputNoiseDoF = oDofs, MaringalLogLikelihood = objs)
ggplot(df, aes(x = OutputNoiseDoF, y = MaringalLogLikelihood)) + geom_line(color = "red") + geom_smooth()
library(urca)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

# Check for cointegration between all explanatory variables
# Johansen test

# Only test unemployment, labour comp, M3 growth, bank rate, interbank rate, and gilt yields
nsIndexes <- c(4, 9, 11, 12, 14, 15)

# 6 time series
#for (m in 2) {
  #combs <- combn(nsIndexes, m)
  #for (cols in 1:ncol(combs)) {
    #j <- combs[, cols]
    #x <- data.frame(inflation.variates.data[, j])
    #print(summary(ca.jo(x, type="trace", K=4, ecdet="none", spec="longrun")))
  #}
#}

x <- data.frame(inflation.variates.data[, c(4, 9, 11, 12, 14, 15)])
print(summary(ca.jo(x, type="trace", K=4, ecdet="none", spec="longrun")))

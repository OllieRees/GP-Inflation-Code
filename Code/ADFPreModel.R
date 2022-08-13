library(tseries)
library(forecast)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

# Check for sationarity
for (i in 2:ncol(inflation.variates.data)) {
  print(colnames(inflation.variates.data)[i])
  data <- inflation.variates.data[, i]
  print(auto.arima(data)$arma[6])
  for (lag in c(1, 4, 12, 25)) print(tseries::adf.test(data, k = lag))
  print(tseries::kpss.test(data))
}
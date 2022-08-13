library(arfima)
library(fracdiff) # for fracdiff and fdSperio/fdGPH
library(forecast) #auto.arima
library(dplyr)

# Data frame cleaning
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)
inflation.basket.data <- read.csv("../Datasets/CPIH_Basket_Quarterly_Reduced.csv", header = T)
rsi.data <- read.csv("../Datasets/RSI.csv", header = T)
trade.data <-  read.csv("../Datasets/TradeValue.csv", header = T)

inflation.variates.data$NumericTime <- 7988:8087/4

# Splitting into training and testing
train <- inflation.variates.data[1:80, ]
test <- inflation.variates.data[81:100, ]

# Get differences and arma orders
d.GPH <- fdGPH(train$CPIH)
train.CPIH.GPH <- diffseries(train$CPIH, d.GPH$d)

d.Sperio <- fdSperio(train$CPIH)
train.CPIH.Sperio <- diffseries(train$CPIH, d.Sperio$d)

train.CPIH.GPH.fit <- auto.arima(train.CPIH.GPH, stepwise = F, approximation = F, parallel = T, num.cores = 20, 
                                 xreg = as.matrix(train[, -c(1, 2, 20)]))
train.CPIH.Sperio.fit <- auto.arima(train.CPIH.Sperio, stepwise = F, approximation = F, parallel = T, num.cores = 20, 
                                    xreg = as.matrix(train[, -c(1, 2, 20)]))

# Get optimal model
fit <- train.CPIH.GPH.fit 
d <- d.GPH
if (train.CPIH.Sperio.fit$aic < train.CPIH.GPH.fit$aic) {
  d <- d.Sperio
  fit <- train.CPIH.Sperio.fit 
}

# Print
print(d)
print(summary(fit))

# Plots
# Plot forecast on unforseen data 
test.forecast <- forecast(fit, h = 20, xreg = as.matrix(test[, -c(1, 2, 20)]))
plot(test.forecast)
print(test.forecast)

# Plot residuals
print(length(test$NumericTime))
plot(x = test$NumericTime, y = sqrt((test$CPIH - test.forecast$mean)^2), xlab = "Time", ylab = "MSE", main = "MSE when Forecasting on Testing Data",
     col = "red", pch = 10)
print(paste("Residual mean: ", mean(sqrt((test$CPIH - test.forecast$mean)^2))))

# Regression errors is the ARMA term and ARIMA error is the noise 
cbind("Regression Errors" = residuals(fit, type="regression"), 
      "ARIMA errors" = residuals(fit, type="innovation")) %>%
  autoplot(facets=TRUE)

checkresiduals(fit)
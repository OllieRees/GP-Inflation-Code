library(arfima)
library(caret)
library(fracdiff) # for fracdiff and fdSperio/fdGPH
library(forecast) #auto.arima
library(ggplot2)

# Data frame cleaning
inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)
inflation.basket.data <- read.csv("../Datasets/CPIH_Basket_Quarterly_Reduced.csv", header = T)

uni.GPH.aic <- c()
uni.Sperio.aic <- c()
multi.GPH.aic <- c()
multi.Sperio.aic <- c()

for (i in 1:100) {
  # Splitting into training and testing
  trainPoints <- createDataPartition(inflation.variates.data$CPIH, p = 0.8, list = F)
  train <- inflation.variates.data[trainPoints, ]
  test <- inflation.variates.data[-trainPoints, ]

  # Get differences and arma orders
  d.GPH <- fdGPH(train$CPIH)
  train.CPIH.GPH <- diffseries(train$CPIH, d.GPH$d)

  d.Sperio <- fdSperio(train$CPIH)
  train.CPIH.Sperio <- diffseries(train$CPIH, d.Sperio$d)

  train.CPIH.GPH.fit <- auto.arima(train.CPIH.GPH, stepwise = F, approximation = F, parallel = T, num.cores = 20)
  train.CPIH.Sperio.fit <- auto.arima(train.CPIH.Sperio, stepwise = F, approximation = F, parallel = T, num.cores = 20)
  
  uni.GPH.aic <- append(uni.GPH.aic, train.CPIH.GPH.fit$aic)
  uni.Sperio.aic <- append(uni.Sperio.aic, train.CPIH.Sperio.fit$aic)
  
  train.CPIH.GPH.fit <- auto.arima(train.CPIH.GPH, stepwise = F, approximation = F, parallel = T, num.cores = 20, 
                                 xreg = as.matrix(train[, 3:19]))
  train.CPIH.Sperio.fit <- auto.arima(train.CPIH.Sperio, stepwise = F, approximation = F, parallel = T, num.cores = 20, 
                                    xreg = as.matrix(train[, 3:19]))

  multi.GPH.aic <- append(multi.GPH.aic, train.CPIH.GPH.fit$aic)
  multi.Sperio.aic <- append(multi.Sperio.aic, train.CPIH.Sperio.fit$aic)
}

print(mean(uni.GPH.aic))
print(mean(uni.Sperio.aic))
print(mean(multi.GPH.aic))
print(mean(multi.Sperio.aic))

means.df <- data.frame(Models = c("Univariate (GPH)", "Univariate (Sperio)", "Multivariate (GPH)", "Multivariate(Sperio)"), 
                       AIC = c(mean(uni.GPH.aic), mean(uni.Sperio.aic), mean(multi.GPH.aic), mean(multi.Sperio.aic)))

ggplot(means.df, aes(x=Models, fill=as.factor(Models), y = AIC)) +  
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "violetred1", "royalblue4", "turquoise1") ) +
  theme(legend.position="none") + 
  ggtitle("The Mean average AIC of various ARFIMA models") + theme(plot.title = element_text(hjust = 0.5))
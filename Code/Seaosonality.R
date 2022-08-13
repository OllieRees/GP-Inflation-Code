library(forecast) #slt

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

nsa.indexes <- c(5, 10, 14, 15, 16, 17, 18)
headers <- c("GDP.Growth", "M1 Growth", "90 day interbank rates", "10 year government backed bonds and securities", "£ effective exchange rate growth",
             "£ effective exchange rate growth by labour cost", "£ effective exchange rate growth by manufacturing CPI")

y <- NULL
for (i in 1:7) {
  j <- nsa.indexes[i]
  x <- ts(inflation.variates.data[, j], frequency = 4)
  season <- stl(x, s.window = 4)
  plot(season, main = paste("Seasonality Plot of ", headers[i]))
  y <- cbind(y, season$time.series[,2] + season$time.series[,3])
}


y <- data.frame(y)
colnames(y) <- headers
z <- apply((abs(inflation.variates.data[, nsa.indexes] - y)), 2, mean)
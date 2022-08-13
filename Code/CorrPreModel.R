library(GGally)
library(RColorBrewer)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)
ggcorr(inflation.variates.data[, 2:19], method = c("everything", "kendall"), label = T, hjust = 1, size = 2)
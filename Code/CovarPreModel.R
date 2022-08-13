library(lattice)
library(RColorBrewer)
library(latticeExtra)

inflation.variates.data <- read.csv("../Datasets/CPIH_Quarterly_Reduced.csv", header = T)

# Covariance matrix
inflation.variates.cov <- cov(inflation.variates.data[, 2:19], method = "kendall")/1000
levelplot(inflation.variates.cov, col.regions = rainbow, main = "Covariance Heatmap of Non-Indexed Independent Vars.",
          scales = list(y = list(at = 1:18, labels = colnames(inflation.variates.data[, 2:19])), 
                        x = list(at = 1:18, labels = 1:18)),
          ylab = "Variables", xlab = "Variable index (bottom to top)"
)

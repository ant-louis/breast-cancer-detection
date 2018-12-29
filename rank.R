# Mahalanobis depth function,
Compute_maha_depth <- function(x, mu, sigma, robust=FALSE){
    if(robust == TRUE){
        #Robust estimation of the covariance matrix
        library(MASS)
        estimator <- cov.rob(x, method = "mcd")
        return(1/(1+ mahalanobis(x, estimator$center, estimator$cov)))
    }
    else{
        1/(1 + mahalanobis(x, mu,sigma))
    } 
}



setwd("~/Documents/INGE/MASTER/1ère\ MASTER/1er\ QUADRI/HDDA/Projects/Breast-cancer-supervised-classification/")

#*******************QUESTION1***************************

# Load quantitative data
data <- read.table("data.csv", header=TRUE, sep=',')
data <- data[,1:9]
attach(data)

# Compute depth
m <- apply(data, 2, mean)
S <- cov(data)
maha_d <- Compute_maha_depth(data,m,S)

# Plot
plot(maha_d)






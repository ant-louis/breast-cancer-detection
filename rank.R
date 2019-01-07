# Set current directory
setwd("~/Documents/INGE/MASTER/1ère\ MASTER/1er\ QUADRI/HDDA/Projects/Breast-cancer-supervised-classification/")


# Mahalanobis depth function
Compute_maha_depth <- function(x, mu=colMeans(x), sigma=cov(x), robust=FALSE){
    if(robust == TRUE){
        #Robust estimation of the covariance matrix
        library(MASS)
        estimator <- cov.rob(x, method = "mcd")
        mu <- estimator$center 
        sigma <- estimator$cov
    }
    1/(1 + mahalanobis(x, mu,sigma))
}

# Load quantitative data
data <- read.table("data.csv", header=TRUE, sep=',')
data <- data[,1:9]
attach(data)


#*******************QUESTION_1***************************
# Compute depth
maha_d <- Compute_maha_depth(data, robust=TRUE)

# Compute PCA
pca <- princomp(data, cor=TRUE)

# Represent data on first principal plane with different ranks
plot(pca$scores[,1], pca$scores[,2],
					main="Data on first principal plane with ranks",
					xlab="First component",
					ylab= "Second component",
					col=topo.colors(length(unique(maha_d)))[as.factor(maha_d)])



#*******************QUESTION_2***************************







#*******************QUESTION_3***************************





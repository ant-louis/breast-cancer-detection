setwd("/home/tom/Documents/Uliège/HDDA/Breast-cancer-supervised-classification")
#setwd("~/Documents/INGE/MASTER/1ère\ MASTER/1er\ QUADRI/HDDA/Projects/Breast-cancer-data-analysis/")

# Data loading
data <- read.table("data.csv", header=TRUE, sep=',')
#Correct data types and values
data$Classification[data$Classification == 1] <- "Healthy"
data$Classification[data$Classification == 2] <- "Cancerous"
data$Classification <-as.factor(data$Classification)
sapply(data[, 1:9], as.numeric)
attach(data)
head(data)

#Check data types
str(data)

#Create training and testing splits
library(caTools)
set.seed(111)
sample = sample.split(data, SplitRatio = .90)
train_data = subset(data, sample == TRUE)
test_data  = subset(data, sample == FALSE)


GLM_complete <- glm(Classification ~ .,
                    data = train_data,
                    control = list(maxit = 50),
                    family=binomial(logit))
summary(GLM_complete)

library(MASS)
# Stepwise computation of the AIC to find the best features
bestAIC <- stepAIC(GLM_complete, direction="both")
print(bestAIC)

# Cross validation
library(glmnet)
CV = cv.glmnet(x=as.matrix(train_data[,-10]),
                y= train_data[,10],
                family="binomial",
                type.measure="class",
                alpha=1,
                nlambda=100)

#Mean and std
print(CV$cvm[which.min(CV$lambda.min)])
print(CV$cvsd[which.min(CV$lambda.min)])

#Coefficients of best model
coef(CV$glmnet.fit, s = CV$lambda.min)

#Plotting residuals and others
library(plotmo) # for plotres

pdf(file="Figures/Coefficients.pdf",title="")
plotres(CV$glmnet.fit,which=1)
dev.off()
pdf(file="Figures/Cumulative_res.pdf",title="")
plotres(CV$glmnet.fit,which=2)
dev.off()
pdf(file="Figures/Residuals.pdf",title="")
plotres(CV$glmnet.fit,which=3)
dev.off()
pdf(file="Figures/QQPlot.pdf",title="")
plotres(CV$glmnet.fit,which=4)
dev.off()


#Predict
pred = predict(CV$glmnet.fit, as.matrix(test_data[,-10]),s=0.001,type="response")
pred_data <- data.frame(proba= pred[,1],class=test_data[,10])
pred_data <- pred_data[order(pred_data$proba, decreasing=FALSE),]
pred_data$rank <- 1:nrow(pred_data)
 
#Plotting prediction vs actual value
library(ggplot2)
library(cowplot)
ggplot(data=pred_data, aes(x=rank, y=proba)) +
  geom_point(aes(color=test_data$Classification), alpha=1, shape=4, stroke=2) +
  xlab("Index") +
  ylab("Predicted probability of getting breast cancer")
 
ggsave("Figures/PredVsClass.pdf")
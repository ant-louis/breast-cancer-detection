L2_reg <- function(S, lambda){
  eig = eigen(S)
  values <- eig$values
  vectors <- eig$vectors
  
  theta_i <- (-values + sqrt(values^2+ 8*lambda))/(4*lambda)
  
  theta_pml <- 0
  theta_pml <- (theta_i * vectors) %*% t(vectors)
  
  res <- list(icov = theta_pml, loglik = log(det(theta_pml)) -sum(diag(S %*% theta_pml)) )
  return(res)
}



setwd("/home/tom/Documents/Uliège/HDDA/Breast-cancer-supervised-classification")
#setwd("~/Documents/INGE/MASTER/1ère\ MASTER/1er\ QUADRI/HDDA/Projects/Breast-cancer-data-analysis/")

# Data loading
data <- read.table("data.csv", header=TRUE, sep=',')
#Correct data types and values
data$Classification[data$Classification == 1] <- "Healthy"
data$Classification[data$Classification == 2] <- "Cancerous"
data$Classification <-as.factor(data$Classification)
attach(data)
View(data)
head(data)
#Check data types
str(data)



GLM_complete <- glm(Classification ~ 
                    data = data,
                    control = list(maxit = 50),
                    family=binomial(logit))
summary(GLM_complete)

#We do a GLM with only the statistically relevant variables (P value < 0.05)
GLM_partial <- glm(Classification ~ BMI + Glucose + Resistin,
                    data = data,
                    control = list(maxit = 50),
                    family=binomial(logit))
summary(GLM_partial)

# --------------------------------------#
# Comparison of the deviance between ----
# two embedded GLM models            ----
# --------------------------------------#
GLM_null <- glm(Classification ~ 1, control = list(maxit = 50), family=binomial(logit))
summary(GLM_null)

anova(GLM_null, GLM_complete, test="Chisq")


# --------------------------------------#
# Selection of variables by AIC     ----
# --------------------------------------#
library(MASS)
# Backward stepwise
stepAIC(GLM_complete)
# Forward stepwise (GLM_null is the NULL model)
stepAIC(GLM_null, scope = Classification ~ Age + BMI + Glucose + Insulin + HOMA 
                    + Leptin + Adiponectin + Resistin + MCP.1)

print(GLM_complete$fitted)
# --------------------------------------#
# Representation of fitted values   ----
# --------------------------------------#
# Representation of the fitted values wrt linear predictor
plot(GLM_complete$linear.predictor, GLM_complete$fitted, col=Classification, pch=16)
legend(x="topleft",col=Classification, pch=16, levels(Classification))

# Representation of the residuals wrt to the fitted values
plot(GLM_complete$fitted, GLM_complete$residuals, col=Classification, pch=16)



predicted.data <- data.frame(
  probability.cancer=GLM_complete$fitted.values,
  class=data$Classification)
 
predicted.data <- predicted.data[
  order(predicted.data$probability.cancer, decreasing=FALSE),]

predicted.data$rank <- 1:nrow(predicted.data)
 
## Lastly, we can plot the predicted probabilities for each sample having
## heart disease and color by whether or not they actually had heart disease
library(ggplot2)
library(cowplot)
ggplot(data=predicted.data, aes(x=rank, y=probability.cancer)) +
  geom_point(aes(color=Classification), alpha=1, shape=4, stroke=2) +
  xlab("Index") +
  ylab("Predicted probability of getting breast cancer")
 
ggsave("heart_disease_probabilities.pdf")












# --------------------------------------#
# Representation of residuals       ----
# --------------------------------------#
plot(GLM_complete$residuals, type='h')
plot(residuals(GLM_complete, type="pearson"), type='h')


# --------------------------------------#
# Prediction for logistic regression ----
# --------------------------------------#
# classification wrt GLM1$fitted  and 0.5
GLM_complete$fitted >= 0.5
table(GLM_complete$fitted >= 0.5, Classification)


# ----------------------------


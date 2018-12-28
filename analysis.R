
setwd("/home/tom/Documents/Uliège/HDDA/Breast-cancer-supervised-classification")
#setwd("~/Documents/INGE/MASTER/1ère\ MASTER/1er\ QUADRI/HDDA/Projects/Breast-cancer-data-analysis/")

# Data loading
data <- read.table("data.csv", header=TRUE, sep=',')
attach(data)

# Data viewer
View(data)


#-----------------------------------------------------------------------------
# Question 2
#-----------------------------------------------------------------------------

# Statistical summaries of the variables
summary(data)

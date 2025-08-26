#####################################################################
## This code demonstrates PCA and SOM in R
## 
## Rachael Blair 
## Created: February 26, 2014
## Edited: March 3, 2015
###################################################################

rm(list = ls())
graphics.off()

library(kohonen)

data(wines)

x11()
boxplot(wines)

# principal components
pc_ex1 <- princomp(wines)
pc_ex2 <- prcomp(wines, center = TRUE, scale = TRUE)
pc_ex3 <- prcomp(scale(wines), center = FALSE, scale = FALSE)

x11()
plot(pc_ex1)

x11()
plot(pc_ex2)

x11()
plot(pc_ex3)

summary(pc_ex3)

# Find the eigenvalues "by hand"
Eigenvalues <- eigen(cov(scale(wines)))
Eigenvalues$values

(pc_ex3$sdev)^2 # sanity check.

quartz()
biplot(pc_ex3)

write.table(wines, file = "wines.txt", sep = "\t", row.names = paste("w",1:177), col.names = colnames(wines))

















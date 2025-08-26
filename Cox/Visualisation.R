library(ggplot2)
#Load data
Boxplot_data <- read_excel("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Results/Boxplot data.xlsx", 
 sheet = "p=200 T1")
#2*2 Multiple graphs
par(mfrow = c(2, 2))
###boxplot for T1 p=200 
### HDnetCox, Lasso, Ridge, ElastcNet
boxplot_temp <- Boxplot_data[,c(9,5,1)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))
boxplot_temp <- Boxplot_data[,c(10,6,2)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))
boxplot_temp <- Boxplot_data[,c(11,7,3)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))
boxplot_temp <- Boxplot_data[,c(12,8,4)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))


#T2 
Boxplot_data <- read_excel("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Results/Boxplot data.xlsx", 
                                sheet = "p=200 T2")
#2*2 Multiple graphs
par(mfrow = c(2, 2))
###boxplot for T2 p=200 
### HDnetCox, Lasso, Ridge, ElastcNet
boxplot_temp <- Boxplot_data[,c(9,5,1)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))
boxplot_temp <- Boxplot_data[,c(10,6,2)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))
boxplot_temp <- Boxplot_data[,c(11,7,3)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))
boxplot_temp <- Boxplot_data[,c(12,8,4)]; colnames(boxplot_temp) <- c("R0.10", "R0.25","R0.40")
boxplot(boxplot_temp, ylim=c(0,100),col=(c("pink","lightblue","lightgreen")))

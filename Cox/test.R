library(survival)
library(SurvMetrics)
time <- rexp(50)
status <- sample(c(0, 1), 50, replace = TRUE)
pre_sp <- runif(50)
t_star <- runif(1)
Brier(Surv(time, status), pre_sp,0)

library(survival)
time <- c(1, 1, 2, 2, 2, 2, 2, 2)
status <- c(0, 1, 1, 0, 1, 1, 0, 1)
predicted <- c(2, 3, 3, 3, 4, 2, 4, 3)
Cindex(Surv(time, status), predicted)
#####


library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

data(veteran)
X_temp <- cbind(time_raw,event_raw,X[,1:10])
#cox <- coxph(Surv(time_raw, event) ~ ., data = X_temp)

cox <- coxph(Surv(time, status) ~ trt + karno +diagtime, data = veteran)
cindex_calcu(time2=veteran$time, status2=veteran$status, X2=veteran[,c(1,5,6)], B2=cox$coefficients)



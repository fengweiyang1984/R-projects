ITERS <- 10
INITIAL_BETA <- 0

GetRiskSet <- function(time_of_interest,
                       entry_vector,
                       time_vector,
                       event_vector) {
  return(which((time_of_interest >= entry_vector) & ((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  #risk set R(T_j) essentially means that the patient hasn’t failed yet or that their censoring date hasn’t passed yet.
  }

#1. X is a n × p matrix representing p covariates for n patients
#2. x_i is a 1 × p vector representing p covariates for patient j
#3. β is a p × 1 vector representing the coefficients for the p covariates
#Calculating the gradient of the likelihood
CoxGradient <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  p <- ncol(Xs)
  
  gradient <- apply(cbind(Ts, Xs), 1, 
                    function(df){
                      
                      df <- matrix(df, nrow = 1)
                      ts <- df[, 1]
                      xs <- df[, 2:ncol(df)]
                      X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event),] %>% 
                        matrix(ncol = ncol(Xs))
                      
                      t1 <- t(X_risk_set) %*% exp(X_risk_set %*% beta)
                      t2 <- sum(exp(X_risk_set %*% beta))
                      
                      return(xs - t1 / t2)
                      
                    }) %>% 
    matrix(nrow = p) %>%
    rowSums()
  return(gradient)
}

#Calculating the variance of the likelihood
CoxVariance <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  
  p <- ncol(Xs)
  
  variance <- apply(cbind(Ts, Xs), 1, 
                    function(df){
                      
                      df <- matrix(df, nrow = 1)
                      ts <- df[, 1]
                      xs <- df[, 2:ncol(df)]
                      X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event),] %>% 
                        matrix(ncol = ncol(Xs))
                      
                      sum1 <- sum(exp(X_risk_set %*% beta))
                      sum2 <- rowSums(t(X_risk_set^2) %*% exp(X_risk_set %*% beta))
                      sum3 <- rowSums(t(X_risk_set) %*% exp(X_risk_set %*% beta))^2
                      
                      return(- (sum1 * sum2 - sum3) / sum1^2)
                      
                    }) %>%
    matrix(nrow = p) %>%
    rowSums()
  
  return(-1 / variance)
  
}

GradientPlot <- function(gradients){
  
  gradient_track <- reshape2::melt(gradients) %>%
    `names<-`(c("iteration", "beta", "gradient"))
  
  p <- ggplot(gradient_track, aes(y = gradient, x = iteration)) +
    geom_line() +
    facet_wrap(~factor(beta))
  
  return(p)
  
}

BetaPlot <- function(store_betas, store_variance, coxph_output){
  
  variance_track <- reshape2::melt(store_variance) %>%
    `names<-`(c("iteration", "variable", "coxph_AS_variance"))
  
  plot_df <- reshape2::melt(store_betas) %>%
    `names<-`(c("iteration", "variable", "coxph_AS_beta")) %>%
    left_join(coxph_output, by = "variable") %>%
    left_join(variance_track, by = c("variable", "iteration")) %>%
    mutate(coxph_AS_lci = coxph_AS_beta - 1.96 * sqrt(coxph_AS_variance),
           coxph_AS_uci = coxph_AS_beta + 1.96 * sqrt(coxph_AS_variance))
  
  ggplot(plot_df) +
    geom_line(aes(y = coxph_AS_beta, x = iteration), size = .8) + 
    geom_line(aes(y = coxph_AS_lci, x = iteration), linetype = "dashed", size = .8) +
    geom_line(aes(y = coxph_AS_uci, x = iteration), linetype = "dashed", size = .8) +
    geom_line(aes(y = coxph_beta, x = iteration, color = "red")) +
    geom_line(aes(y = coxph_lci, x = iteration), linetype = "dashed", color = "red") +
    geom_line(aes(y = coxph_uci, x = iteration), linetype = "dashed", color = "red") +
    facet_wrap(~factor(variable_name)) +
    ylab("beta estimate") +
    theme(legend.position = "none")
  
}

coxph_AS <- function(formula, data){
  
  # Fit using coxph()
  model <- coxph(formula, data)
  cox_beta <- data.frame(variable = 1:length(coef(model)), 
                         coxph_beta = coef(model))
  cox_ci <- cbind(1:length(coef(model)), confint(model)) %>% 
    `colnames<-`(c("variable", "coxph_lci", "coxph_uci")) %>%
    data.frame()
  coxph_output <- left_join(cox_beta, cox_ci, by = "variable") %>%
    mutate(variable_name = names(coef(model)))
  
  # Load data
  x <- data %>% select(all.vars(formula[-1])) %>% as.matrix
  lhs <- setdiff(all.vars(formula), all.vars(formula[-1]))
  
  if (length(lhs) == 1) { # no censoring
    times <- data %>% select(all.vars(formula)[1]) %>% as.matrix
    event <- rep(1, dim(times)[1]) %>% as.matrix
    entry <- rep(0, dim(times)[1]) %>% as.matrix
  } else if (length(lhs) == 2) { # censoring but no entry time specified
    times <- data %>% select(all.vars(formula)[1]) %>% as.matrix
    event <- data %>% select(all.vars(formula)[2]) %>% as.matrix
    entry <- rep(0, dim(times)[1]) %>% as.matrix
  } else if (length(lhs) == 3) { # censoring but no entry time specified
    times <- data %>% select(all.vars(formula)[2]) %>% as.matrix
    event <- data %>% select(all.vars(formula)[3]) %>% as.matrix
    entry <- data %>% select(all.vars(formula)[1]) %>% as.matrix
  }
  # Initialize matrices
  store_gradient <- matrix(NA, nrow = ITERS, ncol = ncol(x))
  store_betas <- matrix(NA, nrow = ITERS, ncol = ncol(x))
  store_variance <- matrix(NA, nrow = ITERS, ncol = ncol(x))
  beta <- matrix(rep(INITIAL_BETA, ncol(x)), nrow = ncol(x))
  
  # Newton-Raphson iterations
  for(i in 1:ITERS){
    store_gradient[i,] <- CoxGradient(beta, x, entry, times, event)
    store_variance[i,] <- CoxVariance(beta, x, entry, times, event)
    store_betas[i,] <- beta <- beta + store_gradient[i,] * store_variance[i,]
  } 
  
  # Plot
  beta_plot <- BetaPlot(store_betas, store_variance, coxph_output)
  gradient_plot <- GradientPlot(store_gradient)
  
  # Final output data.frame
  coxph_output$coxph_AS_beta <- store_betas[i,]
  coxph_output$coxph_AS_lci <- store_betas[i,] - 1.96 * sqrt(store_variance[i,])
  coxph_output$coxph_AS_uci <- store_betas[i,] + 1.96 * sqrt(store_variance[i,])
  
  coxph_output <- coxph_output %>%
    transmute(variable = variable_name,
              beta = coxph_beta,
              beta_AS = coxph_AS_beta,
              LCI = coxph_lci,
              LCI_AS = coxph_AS_lci,
              UCI = coxph_uci,
              UCI_AS = coxph_AS_uci)
  
  return(list(beta_plot = beta_plot,
              gradient_plot = gradient_plot,
              model_output = coxph_output))
}
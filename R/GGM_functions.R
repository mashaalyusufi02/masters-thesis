library(mombf)
library(mvtnorm)
library(tibble)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(ggplot2)
library(pROC)
library(metRology)
library(MASS)      
library(sn)

### df MCMC
run_ggm_analysis_t <- function(p,n,df,w_slab, s_slab = 1, warmup = 5000, iterations = 10000) {
  
  Omega_0 <- diag(x=c(rep(1,p)),nrow=p,ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      if (abs(i-j) == 1){
        Omega_0[i,j] <- 0.5
      }}}
  Omega_0 <- as.matrix(Omega_0)
  Sigma_0 <- solve(Omega_0)
  
  Y <- rmvnorm(n, rep(0, p), Sigma_0)
  
  F_Y <- matrix(NA, nrow = n, ncol = p)
  
  for (i in 1:p){
    F_Y[,i] <- pnorm(Y[,i],0, sqrt(diag(Sigma_0)[i]))
  }
  
  Y_t <- qt.scaled(F_Y,df,sqrt((df-2)/df))
  
  
  Omega_MLE <- solve(cov(Y_t))
  
  posteriorSample <- modelSelectionGGM(
    Y_t, scale = FALSE,
    niter = iterations + warmup,
    burnin = warmup,
    priorCoef = normalidprior(s_slab^2),
    priorModel = modelbinomprior(w_slab),
    priorDiag = exponentialprior(lambda = 1)
    )
  
  graph_inclusion_probabilities <- matrix(NA, nrow = p, ncol = p)
  graph_inclusion_probabilities[upper.tri(graph_inclusion_probabilities, diag = TRUE)] <- posteriorSample$margpp
  
  #Omega_posterior_sample <- array(NA, dim = c(iterations, p, p))
  #for (i in 1:iterations) {
  #  Omega_posterior_sample[i, , ][upper.tri(Omega_posterior_sample[i, , ], diag = TRUE)] <- posteriorSample$postSample[i, ]
  #}
  
  
  
  
  return(list(n=n,graph_inclusion_probabilities=graph_inclusion_probabilities))
}




#skewed distribution

run_ggm_analysis_skew <- function(p,n,skew_param,w_slab, s_slab = 1, warmup = 2000, iterations = 5000) {
  
  
  Omega_0 <- diag(x=c(rep(1,p)),nrow=p,ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      if (abs(i-j) == 1){
        Omega_0[i,j] <- 0.5
      }}}
  Omega_0 <- as.matrix(Omega_0)
  Sigma_0 <- solve(Omega_0)
  
  Y <- rmvnorm(n, rep(0, p), Sigma_0) 
  
  F_Y <- matrix(NA, nrow = n, ncol = p)
  Y_skew <- matrix(NA, nrow=n,ncol=p)
  
  
  b = sqrt(2/pi)
  delta = skew_param/(1+skew_param^2)
  omega = 1/(sqrt(1-(b*delta)^2))
  xi = -b*delta/(sqrt(1-(b*delta)^2))
  
  for (i in 1:p){
    F_Y[,i] <- pnorm(Y[,i],0, sqrt(diag(Sigma_0)[i]))
    Y_skew[,i] <- qsn(F_Y[,i],xi,omega,skew_param)
  }
  
  
  
  
  
  #Omega_MLE <- solve(cov(Y_skew))
  
  suppressMessages(
    invisible(
      capture.output(
        posteriorSample <- modelSelectionGGM(
          Y_skew, scale = FALSE,
          niter = iterations + warmup,
          burnin = warmup,
          priorCoef = normalidprior(s_slab^2),
          priorModel = modelbinomprior(w_slab),
          priorDiag = exponentialprior(lambda = 1)
        ),
        file = NULL
      )
    )
  )
  
  graph_inclusion_probabilities <- matrix(NA, nrow = p, ncol = p)
  graph_inclusion_probabilities[upper.tri(graph_inclusion_probabilities, diag = TRUE)] <- posteriorSample$margpp
  
#  Omega_posterior_sample <- array(NA, dim = c(iterations, p, p))
#  for (i in 1:iterations) {
#    Omega_posterior_sample[i, , ][upper.tri(Omega_posterior_sample[i, , ], diag = TRUE)] <- posteriorSample$postSample[i, ]
#  }
  
  
  
  
  return(list(n=n,graph_inclusion_probabilities=graph_inclusion_probabilities))
}






compute_TPR_FPR_by_sample_loop <- function(results, threshold, type, Omega_0) {
  
  #tpr_fpr_data contains all scoring metrics for all 9 sample sizes (n) for each seed value
  tpr_fpr_data <- tibble(
    seed = character(),
    n = numeric(),
    TPR = numeric(),
    FPR = numeric(),
    FDR = numeric(),
    type = character(),
    threshold = character(),
    AUC = numeric(),
    LS = numeric(),
    BS = numeric()
  )
  
  #roc_list is created to store the AUC values of each seed and each sample size.
  #It isn't used eventually.
  roc_list <- list()
  
  # Loop over each seed (outer list) and each n(inner list)
  for (seed_name in names(results)) {
    seed_results <- results[[seed_name]]
    
    # Loop over each sample size (20,30,40,50..,100)
    for (i in seq_along(seed_results)) {
      #Omega_0 <- seed_results[[i]]$Omega_0
      probabilities <- seed_results[[i]]$graph_inclusion_probabilities
      index <- upper.tri(Omega_0, diag = FALSE)
      
      # Calculate TP, FP, FN, TN
      TP <- sum(Omega_0[index] > 0 & probabilities[index] >= threshold)
      FN <- sum(Omega_0[index] > 0 & probabilities[index] < threshold)
      #FP <- sum((Omega_0[index] == 0) & (probabilities[index] >= (1 - threshold)))
      #TN <- sum((Omega_0[index] == 0) & (probabilities[index] < (1 - threshold)))
      FP <- sum((Omega_0[index] == 0) & (probabilities[index] >= threshold))
      TN <- sum((Omega_0[index] == 0) & (probabilities[index] < threshold))
      
      
      roc_obj <- roc(
        response = as.numeric(Omega_0[index] > 0),
        predictor = probabilities[index]
      )
      auc_val <- as.numeric(auc(roc_obj))
      
      
      #calculate the Log Score
      predval <- probabilities[index]
      truval <- Omega_0[index] > 0 
      LS <- sum(log(predval[truval])) + sum(log(1 - predval[!truval]))

      BS <- mean((probabilities[index] - ifelse(Omega_0[index] > 0, 1, 0))^2)
      
      # Append each row to the tibble tpr_fpr_data
      tpr_fpr_data <- tpr_fpr_data %>% add_row(
        seed = seed_name,
        n = seed_results[[i]]$n,
        TPR = TP / (TP + FN),
        FPR = FP / (FP + TN),
        FDR = ifelse(FP + TP > 0, FP/(FP+TP), 0),
        #type is student-t df=3, student-t df=5, student-t df=100
        type = type,
        threshold = as.character(threshold),
        AUC = auc_val,
        LS = LS,
        BS = BS
      )
      
 
    }
  }
  
  return(tpr_fpr_data)
}






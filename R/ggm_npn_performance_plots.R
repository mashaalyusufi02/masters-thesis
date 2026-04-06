#setwd("~/Honours Coppula Gaussian Graphical Model")
source("GGM_Gt-distributedmodels_functions.R")
ls() 
library(ggplot2)
library(dplyr)
library(pROC)
library(huge)




results_t100 <- readRDS("GGM_tdist_100dof_GGM.rds")
results_t3 <- readRDS("GGM_tdist_3dofGGM.rds")
results_t2.1 <- readRDS("GGM_tdist_2.1dofGGM.rds")


results_skew_0 <- readRDS("GGM_skewresults_standard_normal_skew_alpha_0.rds")
results_skew_2 <- readRDS("GGM_skew_alpha_pos2_GGM.rds")
results_skew_5 <- readRDS("GGM_skew_skew_alpha_pos5_GGM.rds")



p <- 10
Omega_0 <- diag(x=c(rep(1,p)),nrow=p,ncol=p)
for (i in 1:p){
  for (j in 1:p){
    if (abs(i-j) == 1){
      Omega_0[i,j] <- 0.5
    }}}


###function calls to generate marginal posterior probabilities for student t distributions with different degress of freedom

Gaussian_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results, threshold = 0.9,"Gaussian", Omega_0)
t100_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results_t100, threshold = 0.9,"Student-t df=100", Omega_0)
t5_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results_t5, 0.9,"Student-t df=5", Omega_0)
t3_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results_t3, 0.9,"Student-t df=3", Omega_0)
t2.1_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results_t2.1, 0.9,"Student-t df=2.1", Omega_0)

###function calls to generate marginal posterior probabilities for skewed distributions with different degress of freedom
skew0_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results_skew_0, 0.9,"Skew Normal alpha=0", Omega_0)
skew2_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results_skew_2, 0.9,"Skew Normal alpha=2", Omega_0)
skew5_TPR_FPR_0.9_all_seeds <- compute_TPR_FPR_by_sample_loop(results_skew_5, 0.9,"Skew Normal alpha=5", Omega_0)

##bind rows of all tibbles for different models into a single tibble 
final_tibble <- bind_rows(t100_TPR_FPR_0.9_all_seeds,t3_TPR_FPR_0.9_all_seeds,t2.1_TPR_FPR_0.9_all_seeds,skew0_TPR_FPR_0.9_all_seeds,skew2_TPR_FPR_0.9_all_seeds,skew5_TPR_FPR_0.9_all_seeds)#, skew_TPR_FPR_0.95_all_seeds)

final_tibble <- final_tibble %>% mutate(type= as.factor(type))

### FPR graphs


means_tibble <- final_tibble %>% group_by(n, type) %>% summarise(
  meanFPR = mean(FPR),
  meanFDR = mean(FDR),
  meanTPR = mean(TPR),
  meanAUC = mean(AUC),
  meanLS = mean(LS),
  meanBS = mean(BS)
)

means_tibble <- means_tibble %>% filter(n != 10)
### FPR graph for misspecified t-dist cases
FPR_t_graph <- means_tibble %>% 
  filter(type %in% c("Student-t df=100","Student-t df=3","Student-t df=2.1"))%>% 
  ggplot(aes(x = n, y = meanFPR, color = type)) +
  geom_line() +
  labs(
    title = "Gaussian Graphical Model False Positive Rate Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "False Positive Rate (FPR)",
    color = "Type"
  ) 
FPR_t_graph <- FPR_t_graph + ylim(0,0.2)
ggsave("FPR_t_graph_ggm_version_presentation_feedback.png",FPR_t_graph)

### FPR graph for misspecified skew cases
FPR_skew_graph <- means_tibble %>% 
  filter(type %in% c("Skew Normal alpha=0","Skew Normal alpha=5","Skew Normal alpha=2"))%>% 
  ggplot(aes(x = n, y = meanFPR, color = type)) +
  geom_line() +
  labs(
    title = "Gaussian Graphical Model False Positive Rate Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "False Positive Rate (FPR)",
    color = "Type"
  ) 
FPR_skew_graph <- FPR_skew_graph + ylim(0,0.2)
ggsave("FPR_skew_graph_ggm_version_presentation_feedback.png",FPR_skew_graph)

###TPR graphs for misspecified t-dist cases
TPR_t_graph <- means_tibble %>% 
  filter(type %in% 
    c("Student-t df=100","Student-t df=3","Student-t df=2.1"))%>% 
  ggplot(aes(x = n, y = meanTPR, color = type)) +
  geom_line() +
  labs(
    title = "Gaussian Graphical Model True Positive Rate Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "True Positive Rate (TPR)",
    color = "Type"
  ) 
TPR_t_graph <- TPR_t_graph + ylim(0,1)
ggsave("TPR_t_graph_ggm_version_presentation_feedback.png",TPR_t_graph)

###TPR graph for misspecified skew cases
TPR_skew_graph <- means_tibble %>% 
  filter(type %in% c("Skew Normal alpha=0","Skew Normal alpha=5","Skew Normal alpha=2"))%>% 
  ggplot(aes(x = n, y = meanTPR, color = type)) +
  geom_line() +
  labs(
    title = "Gaussian Graphical Model True Positive Rate Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "True Positive Rate (TPR)",
    color = "Type"
  ) 
TPR_skew_graph <- TPR_skew_graph + ylim(0,1)
ggsave("TPR_skew_graph_ggm_version_presentation_feedback.png",TPR_skew_graph)


### FDR graphs
###FDR graphs for misspecified t-dist cases
FDR_t_graph <- means_tibble %>% filter(type %in% c("Student-t df=100","Student-t df=3","Student-t df=2.1"))%>% ggplot(aes(x = n, y = meanFDR, color = type)) +geom_line() + labs(
    title = "Gaussian Graphical Model FDR Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "False Discovery Rate (FDR)",
    color = "Type"
  ) 
FDR_t_graph <- FDR_t_graph + ylim(0,0.5)
ggsave("FDR_t_graph_ggm_version_presentation_feedback.png",FDR_t_graph)

###FDR graph for misspecified skew cases
FDR_skew_graph <- means_tibble %>% filter(type %in% c("Skew Normal alpha=0","Skew Normal alpha=5","Skew Normal alpha=2")) %>% ggplot(aes(x = n, y = meanFDR, color = type)) + geom_line() + labs(
    title = "Gaussian Graphical Model False Discovery Rate Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "False Discovery Rate (FDR)",
    color = "Type"
  ) 
FDR_skew_graph <- FDR_skew_graph + ylim(0,0.5)
ggsave("FDR_skew_graph_ggm_version_presentation_feedback.png",FDR_skew_graph)

### AUC graphs

###AUC graphs for misspecified t-dist cases
AUC_t_graph <- means_tibble %>% 
  filter(type %in% c("Student-t df=100","Student-t df=3","Student-t df=2.1","Gaussian"))%>% 
  ggplot(aes(x = n, y = meanAUC, color = type)) +geom_line() + labs(
    title = "Gaussian Graphical Model AUC Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "Area Under Curve (AUC)",
    color = "Type"
  ) 
AUC_t_graph <- AUC_t_graph + ylim(0,1)
ggsave("AUC_t_graph_ggm_version_presentation_feedback.png",AUC_t_graph)
###AUC graph for misspecified skew cases
AUC_skew_graph <- means_tibble %>% filter(type %in% c("Skew Normal alpha=0","Skew Normal alpha=5","Skew Normal alpha=2")) %>% 
  ggplot(aes(x = n, y = meanAUC, color = type)) + geom_line() + labs(
    title = "Gaussian Graphical Model AUC Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "Area Under Curve (AUC)",
    color = "Type"
  ) 
AUC_skew_graph <- AUC_skew_graph + ylim(0,1)
ggsave("AUC_skew_graph_ggm_version_presentation_feedback.png",AUC_skew_graph)
### LS graphs

###LS graphs for misspecified t-dist cases
LogScore_t_graph <-means_tibble %>% 
  filter(type %in% c("Student-t df=100","Student-t df=3","Student-t df=2.1"))%>% 
  ggplot(aes(x = n, y = meanLS, color = type)) +geom_line() + labs(
    title = "Gaussian Graphical Model LS Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "Log Score",
    color = "Type"
  ) 

ggsave("LogScore_t_graph_ggm_version_presentation_feedback.png",LogScore_t_graph)

###LS graph for misspecified skew cases
LogScore_skew_graph <- means_tibble %>% filter(type %in% c("Skew Normal alpha=0","Skew Normal alpha=5","Skew Normal alpha=2")) %>% 
  ggplot(aes(x = n, y = meanLS, color = type)) + geom_line() + labs(
    title = "Misspecified Skew Margins LS Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "Log Score",
   color = "Type"
  ) 

ggsave("LogScore_skew_graph_ggm_version_presentation_feedback.png",LogScore_skew_graph)


BrierScore_t_graph <-means_tibble %>%
  filter(type %in% c("Student-t df=100","Student-t df=3","Student-t df=2.1"))%>% 
  ggplot(aes(x = n, y = meanBS, color = type)) +geom_line() + labs(
    title = "Gaussian Graphical Model Brier Score Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "Brier Score",
    color = "Type"
  ) 

BrierScore_t_graph + ylim(0,0.09)
ggsave("BrierScore_t_graph_ggm_version_presentation_feedback.png",BrierScore_t_graph)



BrierScore_skew_graph <- means_tibble %>% filter(type %in% c("Skew Normal alpha=0","Skew Normal alpha=5","Skew Normal alpha=2")) %>% 
  ggplot(aes(x = n, y = meanBS, color = type)) + geom_line() + labs(
    title = "Gaussian Graphical Model Brier Score Averaged across 100 samples",
    x = "Sample Size (n)",
    y = "Brier Score",
    color = "Type"
  ) 

BrierScore_skew_graph + ylim(0,0.09)
ggsave("BrierScore_skew_graph_ggm_version_presentation_feedback.png",BrierScore_skew_graph)

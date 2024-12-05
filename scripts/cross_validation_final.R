library(tidyverse)
library(rstan)
library(loo)
library(bayesplot)
library(caret)
library(cmdstanr)
library(lme4)

setwd('~/Dropbox/UNC/Fall2024/BIOS784/Final_Project/BIOS-784-Project/scripts/')

get_folds = function(df, kfolds = 5) {
  n <- nrow(df)
  grps <- unique(df$region_grp)
  ngrps <- length(grps)
  df$fold <- NA
  for (y in c(0,1)) {
    for (g in grps) {
      inds_gy = which(df$region_grp == g & df$resistance == y)
      folds = createFolds(inds_gy, k=kfolds, list=T)
      for (fold_idx in seq_along(folds)) {
        df$fold[inds_gy[folds[[fold_idx]]]] <- fold_idx
      }
    }
  }
  return(df)
}
get_summ_metrics_bin = function(y, yhat, probs=NULL) {
  TP <- sum(yhat[y==1])
  TN <- sum(yhat[y==0]==0)
  FP <- sum(yhat[y==0]==1)
  FN <- sum(yhat[y==1]==0)
  n <- length(y)
  accuracy <- (TP + TN) / n
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  precision <- TP / (TP + FP)
  if (TP + FP + FN > 0) {
    f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  } else {
    f1_score <- NA
  }
  RMSE <- NA
  if (!is.null(probs)) {
    RMSE <- sqrt(mean((y-probs)^2))
  }
  res <- data.frame(accuracy = accuracy,
                    F1_score = f1_score,
                    sensitivity = sensitivity,
                    specificity = specificity,
                    RMSE = RMSE)
  return(res)
}
competing_mods_predict = function(df_train, df_test) {
  df_train$region_grp = factor(df_train$region_grp)
  df_test$region_grp = factor(df_test$region_grp)

  mod_main = glm(resistance ~ Monoclonal + scale(year) + region_grp, data = df_train, family = stats::binomial)
  mod_int = glm(resistance ~ Monoclonal + scale(year) + Monoclonal*region_grp, data = df_train, family = stats::binomial)
  mod_ri = glmer(resistance ~ Monoclonal + scale(year) + (1|region_grp), data = df_train, family = stats::binomial)
  mod_rs = glmer(resistance ~ Monoclonal + scale(year) + (Monoclonal+scale(year)|region_grp), data = df_train, family = stats::binomial)

  probs_main = predict(mod_main, newdata = df_test, type = "response")
  probs_int = predict(mod_int, newdata = df_test, type = "response")
  probs_ri = predict(mod_ri, newdata = df_test, type = "response")
  probs_rs = predict(mod_rs, newdata = df_test, type = "response")
  
  yhat_main = ifelse(probs_main >= .5, 1, 0)
  yhat_int = ifelse(probs_int >= .5, 1, 0)
  yhat_ri = ifelse(probs_ri >= .5, 1, 0)
  yhat_rs = ifelse(probs_rs >= .5, 1, 0)
  
  yobs = df_test$resistance
  res = rbind(data.frame(Model="Fixed Main", get_summ_metrics_bin(yobs, yhat_main, probs_main)),
        data.frame(Model = "Fixed w/ Interaction", get_summ_metrics_bin(yobs, yhat_int, probs_int)),
        data.frame(Model = "Random Intercept", get_summ_metrics_bin(yobs, yhat_ri, probs_ri)),
        data.frame(Model = "Random Slope", get_summ_metrics_bin(yobs, yhat_rs, probs_rs)))
  return(res)
}
cross_validation = function(df, kfolds = 5, nrep=1, seed=123) {
  set.seed(seed)
  res = data.frame()
  for (r in 1:nrep) {
    print(paste0('Iter: ',r))
    df_cv = get_folds(df, kfolds)
    for (k in 1:kfolds) {
      print(paste0('Fold=',k))
      df_train_k = df_cv %>% filter(fold != k)
      df_test_k = df_cv %>% filter(fold == k)
      res_rk = competing_mods_predict(df_train_k, df_test_k)
      res_rk$Fold = k
      res_rk$Iter = r
      res = rbind(res, res_rk)
    }
  }
  return(res)
}
get_outcome_data = function(df, outcome = "CQ", group = "country", n_cutoff=100, sd_y_cutoff = 0) {
  df = df %>% filter(drug == outcome)
  if (group == "country") df$region_grp = as.integer(factor(df$country))
  df = df %>% group_by(region_grp) %>% mutate(ngrp = n(), sd_y = sd(resistance)) %>% filter(ngrp > n_cutoff & sd_y > sd_y_cutoff)
  df$region_grp = as.integer(factor(df$region_grp))
  return(df)
}

df = read.csv('../data/pf7_prediction_orig.csv')
# df = read.csv('../data/pf7_prediction.csv')

kfolds = 10
nrep = 5
seed = 12
cv_art = df %>% get_outcome_data(outcome = "ART") %>% cross_validation(kfolds = kfolds, nrep = nrep, seed = seed)
write.csv(cv_art, '../cv_results/cv_res_art_10foldB.csv')
cv_cq = df %>% get_outcome_data(outcome = "CQ") %>% cross_validation(kfolds = kfolds, nrep = nrep, seed = seed)
write.csv(cv_cq, '../cv_results/cv_res_cq_10foldB.csv')

library(RColorBrewer)
make_plots_cv = function(res) {
  res2 = res %>%
    mutate(Model = factor(Model, levels = c("Fixed w/ Interaction", "Random Slope", "Fixed Main", "Random Intercept")))
    
  p_acc = res2 %>%
    ggplot(aes(x = Model, y = accuracy, fill = Model))+
    geom_boxplot()+
    scale_fill_brewer(palette = "Set1")+
    labs(fill = "Model")+
    theme_bw(base_size=18)+
    ylab('Accuracy')+
    xlab('Model')
  
  p_f1 = res2 %>%
    ggplot(aes(x = Model, y = F1_score, fill = Model))+
    geom_boxplot()+
    scale_fill_brewer(palette = "Set1")+
    labs(fill = "Model")+
    theme_bw(base_size=18)+
    ylab('F1 Score')+
    xlab('Model')
  
  return(list(p_acc = p_acc, p_f1 = p_f1))
}

cv_art = read.csv('../cv_results/cv_res_art_10foldB.csv')
cv_cq = read.csv('../cv_results/cv_res_cq_10foldB.csv')

make_plots_cv(cv_art)$p_acc+ggtitle('ART')+
  theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  ggtitle('ART Resistance (10-fold CV)')
make_plots_cv(cv_art)$p_f1+ggtitle('ART')

make_plots_cv(cv_cq)$p_acc+ggtitle('CQ')+
  theme(legend.position = "none",axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  ggtitle('CQ Resistance (10-fold CV)')
make_plots_cv(cv_cq)$p_f1+ggtitle('CQ')
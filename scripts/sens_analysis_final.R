library(tidyverse)
library(rstan)
library(loo)
library(bayesplot)
library(caret)
library(cmdstanr)
library(lme4)

setwd('~/Dropbox/UNC/Fall2024/BIOS784/Final_Project/BIOS-784-Project/scripts/')

split_data = function(df, p = .2) {
  n <- nrow(df)
  grps <- unique(df$region_grp)
  ngrps <- length(grps)
  train_idx <- c()
  for (g in grps) {
    inds_g <- which(df$region_grp == g)
    df_g <- df[inds_g,]
    train_idx_g <- createDataPartition(df_g$resistance, p=p, list = F)
    train_idx <- c(train_idx, inds_g[train_idx_g])
  }
  train_data <- df[train_idx, ]
  test_data <- df[-train_idx, ]
  return(list(df_train = train_data, df_test = test_data))
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
sens_analysis = function(df, p = .1, nrep=3, seed=123) {
  set.seed(seed)
  res = data.frame()
  for (r  in 1:nrep) {
    print(r)
    df_split = split_data(df, p)
    res_r = competing_mods_predict(df_split$df_train, df_split$df_test)
    res = rbind(res, res_r)
  }
  p_acc = res %>% ggplot(aes(y = accuracy, fill = factor(Model)))+
    geom_boxplot()+
    theme_bw()+
    labs(fill = "Model")+
    ylab("Accuracy")
  p_f1 = res %>% ggplot(aes(y = F1_score, fill = factor(Model)))+
    geom_boxplot()+
    theme_bw()+
    labs(fill = "Model")+
    ylab("F1 Score")
  return(list(res=res, p1=p_acc, p2=p_f1))
}
get_outcome_data = function(df, outcome = "CQ", group = "country", n_cutoff=100, sd_y_cutoff = 0) {
  df = df %>% filter(drug == outcome)
  if (group == "country") df$region_grp = as.integer(factor(df$country))
  df = df %>% group_by(region_grp) %>% mutate(ngrp = n(), sd_y = sd(resistance)) %>% filter(ngrp > n_cutoff & sd_y > sd_y_cutoff)
  df$region_grp = as.integer(factor(df$region_grp))
  return(df)
}
sens_analysis_full = function(df, outcome = "CQ", ps = 1:4/5, nrep = 10, seed = 123, save = T) {
  print(outcome)
  df_o = get_outcome_data(df, outcome)
  print(dim(df_o))
  res = data.frame()
  for (p in ps) {
    print(p)
    res_p = sens_analysis(df_o, p = p, nrep = nrep, seed = seed)
    res_p$res$Train_Size = p
    res = rbind(res, res_p$res)
  }
  if (save == T) {
    write.csv(res, paste0('../sens_analysis_results/res_',tolower(outcome),'B.csv'))
  }
  return(res)
}

df = read.csv('../data/pf7_prediction_orig.csv')

res_art = sens_analysis_full(df, outcome = "ART")
res_cq = sens_analysis_full(df, outcome = "CQ")

library(RColorBrewer)
make_plots = function(res, k=4) {
  p_acc = ggplot(res, aes(x = Train_Size, y = accuracy, col = factor(Model), fill = factor(Model)))+
    geom_smooth(data = res %>% filter(Model %in% c("Random Slope", "Fixed w/ Interaction")), 
                formula = y ~ s(x, k=k), method = "gam", se = T, alpha = .25)+
    geom_smooth(data = res %>% filter(Model %in% c("Random Intercept", "Fixed Main")), 
                formula = y ~ s(x, k=k), method = "gam", se = F, linetype = "dashed")+
    scale_color_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    labs(col = "Model", fill = "Model")+
    theme_bw(base_size=18)+
    ylab('Accuracy')+
    xlab('Training Size')
  
  p_f1 = ggplot(res, aes(x = Train_Size, y = F1_score, col = factor(Model), fill = factor(Model)))+
    geom_smooth(data = res %>% filter(Model %in% c("Random Slope", "Fixed w/ Interaction")), 
                formula = y ~ s(x, k=k), method = "gam", se = T, alpha = .25)+
    geom_smooth(data = res %>% filter(Model %in% c("Random Intercept", "Fixed Main")), 
                formula = y ~ s(x, k=k), method = "gam", se = F, linetype = "dashed")+
    scale_color_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    labs(col = "Model", fill = "Model")+
    theme_bw(base_size=18)+
    ylab('F1 Score')+
    xlab('Training Size')
  
  return(list(p_acc = p_acc, p_f1 = p_f1))
}

res_art = read.csv('../sens_analysis_results/res_artB.csv')
res_cq = read.csv('../sens_analysis_results/res_cqB.csv')

make_plots(res_art)$p_acc+ggtitle('ART Resistance')+theme(legend.position = "none")
make_plots(res_art)$p_f1+ggtitle('ART')

make_plots(res_cq, k = 4)$p_acc+ggtitle('CQ Resistance')+theme(legend.position = "none")
make_plots(res_cq)$p_f1+ggtitle('CQ')
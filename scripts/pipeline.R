library(tidyverse)
library(lme4)
library(mice)
library(caret)
library(logisticPCA)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rlist)

# Data Cleaning ------------------------------------------------------------

#setwd('~/Dropbox/UNC/Fall2024/BIOS784/Final_Project/')
fws_data <- read.csv('../data/pf7_fws.txt', sep = "\t")
fws_data$sample_id <- fws_data$Sample
metadata <- read.csv('../data/Pf7-samples.csv')
df_merge <- inner_join(fws_data, metadata, by = "sample_id") %>% 
  select(-Sample)
df_merge <- df_merge %>%
  mutate(across(ends_with("resistant"), ~ case_when(
    . == "sensitive" ~ 0,
    . == "resistant" ~ 1,
    . == "undetermined" ~ NA_real_
  ))) %>%
  mutate(Monoclonal = ifelse(Fws >= .95, 1, 0))

# include spatial location
df_merge <- df_merge %>%
  geocode(country = country, method = "osm") %>%
  mutate(country = as.factor(country))

write_csv(df_merge, '../data/pf7_merged.csv')

# Data Processing ---------------------------------------------------------

df_model <- read.csv('../data/pf7_merged.csv')

# transform Fws
ep <- .001
df_model <- df_model %>%
  mutate(Fws_sc = Fws*(1-2*ep)+ep) %>%
  mutate(Fws_logit = log(Fws_sc/(1-Fws_sc)))
hist(df_model$Fws)
hist((df_model$Fws_logit))

# impute missing resistances
set.seed(123)
outcome_vars <- grep("resistant$", names(df_model), value = T)
df_toImpute <- df_model %>%
  select(all_of(outcome_vars), Monoclonal, country)
mice_res <- mice(df_toImpute, method = "pmm", m = 5, seed = 123)
df_complete <- complete(mice_res, 1)
df_imputed <- df_model
df_imputed[outcome_vars] <- df_complete[outcome_vars]
write.csv(df_imputed, '../data/pf7_imputed.csv')

# logistic PCA on imputed outcomes
df_pca <- df_imputed
df_resistances <- df_pca[outcome_vars]
log_pca_res <- logisticPCA(df_resistances, k = 2, m = 10)
df_pca$PC1 <- log_pca_res$PCs[,1]
df_pca$PC2 <- log_pca_res$PCs[,2]
write.csv(df_pca, '../data/pf7_pca.csv')

# cluster by country
df_pca %>%
  ggplot(aes(x = PC1, y = PC2, color = country))+
  geom_point(position = position_jitter(width=2,height=2), alpha = .25)+
  theme_minimal()

set.seed(23)
nclusters <- 6
#x_cluster <- df_pca[c("Monoclonal","PC1","PC2","long","lat")]
x_cluster <- df_pca[c(outcome_vars,"Monoclonal","long","lat")]
kmeans_res <- kmeans(x = x_cluster, centers = nclusters)
df_pca$cluster <- kmeans_res$cluster
all_clusters <- sort(unique(df_pca$cluster))
df_clusters <- df_pca %>%
  group_by(country) %>%
  mutate(
    assigned_cluster = names(which.max(prop.table(table(factor(cluster, levels = all_clusters)))))
  ) %>%
  mutate(region_grp = factor(assigned_cluster)) %>%
  select(-cluster)
df_clusters_byCountry <- df_clusters %>%
  group_by(country) %>%
  summarise(cluster = unique(assigned_cluster)[1])

world <- ne_countries(scale = "medium", returnclass = "sf")
map_data <- world %>%
  mutate(country = admin) %>%
  left_join(df_clusters_byCountry, by = "country")
ggplot(data = map_data) +
  geom_sf(aes(fill = cluster)) +
  scale_fill_brewer(palette = "Set3", na.translate = FALSE) + # Use a color palette
  labs(
    title = "Country Clusters",
    fill = "Cluster"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

write.csv(df_clusters, '../data/pf7_clustered.csv')

df_final <- df_clusters
# df_final <- df_clusters %>%
#   mutate(nresistances = CQresistant + MQresistant + ARTresistant + PPQresistant + PYRresistant + SDXresistant) %>%
#   mutate(resistant1 = ifelse(nresistances >= 1, 1, 0)) %>%
#   mutate(resistant2 = ifelse(nresistances >= 2, 1, 0)) %>%
#   mutate(resistant3 = ifelse(nresistances >= 3, 1, 0)) %>%
#   mutate(resistant4 = ifelse(nresistances >= 4, 1, 0)) %>%
#   mutate(resistant5 = ifelse(nresistances >= 5, 1, 0)) %>%
#   mutate(resistant6 = ifelse(nresistances >= 6, 1, 0))
# 
# write.csv(df_final, '../data/pf7_prediction.csv')
  
df_long = df_final %>%
  pivot_longer(cols = ends_with("resistant"), names_to = "drug", values_to = "resistance") %>%
  mutate(resistance_f = factor(ifelse(resistance == 1, "Yes", "No"), levels = c("Yes","No"))) %>%
  #mutate(resistance_f = factor(resistance, levels = c(1,0))) %>%
  mutate(drug = factor(sub("resistant.*$", "", drug))) %>%
  mutate(region_grp = factor(region_grp, levels = 1:6))


# monoclonal frequencies
df_long %>%
  ggplot(aes(x = factor(region_grp), fill = factor(Monoclonal)))+
  geom_bar(position = "fill")+
  ylab("Proportion")+
  xlab("Country")+
  labs(fill = "Resistant")+
  ggtitle('Relative frequency of drug resistance')


# resistance frequencies
df_long %>%
  ggplot(aes(x = drug, fill = resistance_f))+
  geom_bar(position = "fill")+
  facet_wrap(vars(region_grp),
             labeller = as_labeller(function(label) paste("Region", label)))+
  ylab("Proportion")+
  xlab("Drug")+
  labs(fill = "Resistant")+
  ggtitle('Relative frequency of drug resistance')

# resistance on Monoclonal trends
df_long %>%
  mutate(Polyclonal = 1-Monoclonal) %>%
  ggplot(aes(x = Polyclonal, y = resistance, col = factor(region_grp)))+
  stat_smooth(method = "glm", se = F, method.args = "binomial")+
  facet_wrap(vars(drug))+
  labs(col = "Region")+
  theme_bw()+
  ylab('Reistance')+
  ggtitle('Fitted resistance probabilities')

df_long %>%
  group_by(region_grp, drug) %>%
  summarise(sd_x = sd(Monoclonal), sd_y = sd(resistance))

# analyze variation in outcome by group
df_grp = df_long %>%
  group_by(drug, region_grp, Monoclonal) %>%
  summarise(n = n(), sd_y = sd(resistance), .groups = "keep")
df_grp_novar = df_grp %>% filter(sd_y == 0)

write.csv(df_long, '../data/pf7_prediction_orig.csv')

df_long_rm = df_long %>%
  group_by(drug, region_grp, Monoclonal) %>%
  mutate(sd_y = sd(resistance)) %>%
  filter(sd_y != 0)

write.csv(df_long_rm, '../data/pf7_prediction.csv')
  

# Fit Full Models ---------------------------------------------------------
df_long_rm = read.csv('../data/pf7_prediction.csv')

fit_models = function(data, outcome=NULL) {
  if (!is.null(outcome)) {
    data_o = data %>% filter(drug == outcome)
  } else {
    data_o = data
  }
  mod_main_o = glm(resistance ~ Monoclonal+scale(year)+region_grp, data = data_o, family = "binomial")
  mod_int_o = glm(resistance ~ scale(year)+Monoclonal*region_grp, data = data_o, family = "binomial")
  mod_ri_o = glmer(resistance ~ scale(year)+Monoclonal + (1|region_grp), data = data_o, family = "binomial")
  mod_rs_o = glmer(resistance ~ scale(year)+Monoclonal + (Monoclonal+scale(year)|region_grp), data = data_o, family = "binomial")
  Model_list = list("Fixed_Main"=mod_main_o,
                    "Fixed_Int"=mod_int_o,
                    "Mixed_RI"=mod_ri_o,
                    "Mixed_RS"=mod_rs_o)
  return(Model_list)
}

res_art = fit_models(df_long_rm, outcome = "ART")
res_cq = fit_models(df_long_rm, outcome = "CQ")
res_mq = fit_models(df_long_rm, outcome = "MQ")
res_ppq = fit_models(df_long_rm, outcome = "PPQ")
res_pyr = fit_models(df_long_rm, outcome = "PYR")
res_sdx = fit_models(df_long_rm, outcome = "SDX")

# get AICs
unlist(lapply(res_art, AIC))
unlist(lapply(res_cq, AIC))
unlist(lapply(res_mq, AIC))
unlist(lapply(res_ppq, AIC))
unlist(lapply(res_pyr, AIC))
unlist(lapply(res_sdx, AIC))



# Prediction --------------------------------------------------------------
fit_models_cv = function(data, outcome=NULL) {
  if (!is.null(outcome)) {
    data_o = data %>% filter(drug == outcome)
  } else {
    data_o = data
  }
  mod_main_o = glm(resistance ~ Monoclonal+region_grp+year, data = data_o, family = "binomial")
  mod_int_o = glm(resistance ~ Monoclonal*region_grp+year, data = data_o, family = "binomial")
  mod_ri_o = glmer(resistance ~ Monoclonal + year + (1|region_grp), data = data_o, family = "binomial")
  mod_rs_o = glmer(resistance ~ Monoclonal + year + (Monoclonal|region_grp), data = data_o, family = "binomial")
  Model_list = list("Fixed_Main"=mod_main_o,
                    "Fixed_Int"=mod_int_o,
                    "Mixed_RI"=mod_ri_o,
                    "Mixed_RS"=mod_rs_o)
  return(Model_list)
}
get_folds = function(data, groups, nfolds = 10) {
  data$group = do.call(paste, c(data[groups], sep="_"))
  data$fold = NA
  unique_groups = unique(data$group)
  for (g in unique_groups) {
    group_indices = which(data$group == g)
    folds = createFolds(group_indices, k=nfolds, list=T)
    for (fold_idx in seq_along(folds)) {
      data$fold[group_indices[folds[[fold_idx]]]] <- fold_idx
    }
  }
  return(data)
}
get_summ_metrics = function(y, probs) {
  yhat <- ifelse(probs >= .5, 1, 0)
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
  RMSE <- sqrt(mean((y-probs)^2))
  res <- data.frame(accuracy = accuracy,
                    F1_score = f1_score,
                    sensitivity = sensitivity,
                    specificity = specificity,
                    RMSE = RMSE)
  return(res)
}
model_predict = function(mod, newdata) {
  probs = predict(mod, newdata, type = "response")
  y = newdata$resistance
  return(get_summ_metrics(y, probs))
}
cross_validation = function(data, outcome, groups, nfolds = 10, nreps=1, seed = 123) {
  pred_res = vector(mode = "list", 4)
  for (r in 1:nreps) {
    print(paste0('Replicate = ', r))
    set.seed(seed + r)
    df_cv = data %>% filter(drug == outcome) %>% get_folds(groups, nfolds)
    for (k in 1:nfolds) {
      print(paste0('Fold = ',k))
      df_train_k = df_cv %>% filter(fold != k)
      df_test_k = df_cv %>% filter(fold == k)
      mods_k = fit_models_cv(df_train_k)
      preds_k = lapply(mods_k, function(x) cbind(model_predict(x, df_test_k), data.frame("fold"=k, "rep"=r)))
      pred_res = lapply(1:4, function(x) {
        rbind(pred_res[[x]], preds_k[[x]])
      })
    }
  }
  names(pred_res) = names(preds_k)
  return(pred_res)
}

df_long_rm = read.csv('../data/pf7_prediction.csv')
groups = c("resistance", "Monoclonal", "region_grp")
unique(df_long_rm$drug)

# Run 10-fold cross-validation
cv_res_art = cross_validation(df_long_rm, "ART", groups, nfolds = 10, nreps=3)
cv_res_cq = cross_validation(df_long_rm, "CQ", groups, nfolds = 10, nreps=3)
cv_res_ppq = cross_validation(df_long_rm, "PPQ", groups, nfolds = 10, nreps=3)
cv_res_pyr = cross_validation(df_long_rm, "PYR", groups, nfolds = 10, nreps=3)
cv_res_sdx = cross_validation(df_long_rm, "SDX", groups, nfolds = 10, nreps=3)
cv_res_mq = cross_validation(df_long_rm, "MQ", groups, nfolds = 10, nreps=3)

# Save results
list.save(cv_res_art, '../cv_results/ARTPrediction.rdata')
list.save(cv_res_cq, '../cv_results/CQPrediction.rdata')
list.save(cv_res_ppq, '../cv_results/PPQPrediction.rdata')
list.save(cv_res_pyr, '../cv_results/PYRPrediction.rdata')
list.save(cv_res_sdx, '../cv_results/SDXPrediction.rdata')
list.save(cv_res_mq, '../cv_results/MQPrediction.rdata')


# Plot Results ------------------------------------------------------------
library(ggpubr)

#cv_res_cq = list.load('../cv_results/CQPrediction.rdata')
df_cv_res_cq = do.call(rbind, lapply(names(cv_res_cq), function(model_name) {
  df = cv_res_cq[[model_name]]
  df$Model = model_name
  return(df)
}))
df_cv_plot_cq = df_cv_res_cq %>%
  mutate(Model = factor(case_when(Model == "Fixed_Int" ~ "Fixed Effects w/ interaction",
                                  Model == "Fixed_Main" ~ "Fixed Effects",
                                  Model == "Mixed_RI" ~ "Random Intercept",
                                  Model == "Mixed_RS" ~ "Random Slope")))

p = df_cv_plot_cq %>%
  ggplot(aes(x = Model, y = F1_score, fill = Model))+
  geom_boxplot()+
  labs(fill = "Model")+
  theme_bw(base_size=18)+
  ylab('F1 Score')+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab('')
stat_test = compare_means(F1_score ~ Model, data = df_cv_plot_cq, method = "t.test")
p = p+stat_compare_means(
  comparisons = list(c("Fixed Effects w/ interaction", "Random Intercept")),
  method = "t.test",
  label = "p.signif"
)
p




# Stan Modeling -----------------------------------------------------------
library(rstan)
library(loo)
library(bayesplot)
library(caret)
library(cmdstanr)

stan_train_mod = function(df_stan, p = .2, chains = 4, niter = 2000, seed = 42) {
  set.seed(seed)
  n <- nrow(df_stan)
  train_idx <- createDataPartition(df_stan$resistance, p = .2, list = F)
  train_data <- df_stan[train_idx, ]
  test_data <- df_stan[-train_idx, ]
  stan_train_data <- list(
    N = nrow(train_data),
    I = length(unique(train_data$region_grp)),
    region = train_data$region_grp,
    x = train_data$Monoclonal,
    y = train_data$resistance
  )
  stan_mod <- cmdstan_model('../stan/hier_model_v3.stan')
  iter_sampling = iter_warmup = round(niter/2)
  stan_fit <- stan_mod$sample(
    data = stan_train_data,
    chains = chains,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    seed = seed, 
    parallel_chains = chains
  )
  res = list(stan_fit = stan_fit,
             train_data = train_data,
             test_data = test_data)
}
stan_eval_mod = function(stan_fit, test_data) {
  stan_test_data <- list(
    N = nrow(test_data),
    I = length(unique(train_data$region_grp)),
    region = test_data$region_grp,
    x = test_data$Monoclonal
  )
  alpha0 <- stan_fit$draws("alpha0", format = "matrix")
  alpha1 <- stan_fit$draws("alpha1", format = "matrix")
  mu0 <- stan_fit$draws("mu0", format = "matrix")
  mu1 <- stan_fit$draws("mu1", format = "matrix")
  
  # predict
  n_test <- nrow(test_data)
  n_samples <- length(alpha0)
  y_pred_test <- matrix(NA, nrow = n_samples, ncol = n_test)
  y_obs_test <- test_data$resistance
  
  for (i in 1:n_samples) {
    eta <- alpha0[i] + alpha1[i] * test_data$Monoclonal +
      mu0[i, test_data$region_grp] + mu1[i, test_data$region_grp] * test_data$Monoclonal
    y_pred_test[i, ] <- rbinom(n_test, size = 1, prob = plogis(eta))
  }
  
  y_pred_probs_test <- apply(y_pred_test,2,mean)
  y_hat_test <- ifelse(y_pred_probs_test >= .5,1,0)
  
  print(mean(y_hat_test==y_obs_test))
  return(list(yobs = y_obs_test, ypred = y_hat_test, y_pred_probs = y_pred_probs_test, ypred_samps = y_pred_test))
}

df_stan <- df_long_rm %>% filter(drug == "CQ")
train_res = stan_train_mod(df_stan, p = .1, chains = 2, niter = 1500, seed = 123)
test_res = stan_eval_mod(train_res$stan_fit, train_res$test_data)



















set.seed(123)
N = 500
df_stan = df_long_rm %>% filter(drug == "CQ")
nobs = nrow(df_stan)
samps = sample(1:nobs, N)
region = as.integer(df_stan$region_grp[samps])
y = df_stan$resistance[samps]
x = df_stan$Monoclonal[samps]
nregions = length(unique(region))

niter <- 2000

stan_data <- list(
  N = N, I = nregions, region = region, x = x, y = y
)

fit_hier <- sampling(
  stan_model(file = "../stan/hier_model_v2.stan"),
  data = stan_data,
  iter = niter, 
  chains = 2,
  seed = 123
)

y_rep = extract(fit_hier, pars = "y_rep")$y_rep
ppc_stat = bayesplot::ppc_stat(y = y, yrep = y_rep)
bayesplot::ppc_dens_overlay(y = y, yrep = y_rep[1:50, ])


# Prediction (10-fold Cross-Validation) --------------------------------------------------------------
df_final = read.csv('../data/pf7_prediction.csv')
data = df_final %>% na.omit() %>%
  mutate(site_id = factor(site_id)) %>%
  mutate(country = factor(country))

remove_zerovar <- function(data, groups) {
  y_var <- groups[1]
  data$group <- do.call(paste, c(data[groups], sep="_"))
  if (length(groups) == 1) {
    data$group2 <- data$group
    data_rm <- data
  } else {
    data$group2 <- sub("^[^_]*_", "", data$group)
    data_rm <- data %>%
      group_by(group2) %>%
      mutate(sd = sd(!!sym(y_var))) %>%
      filter(sd != 0) %>%
      ungroup()
  }
  return(as.data.frame(data_rm))
}
get_folds <- function(data, groups, nfolds=10, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  data$group <- do.call(paste, c(data[groups], sep="_"))
  data$fold <- NA
  unique_groups <- unique(data$group)
  for (g in unique_groups) {
    group_indices <- which(data$group == g)
    folds <- createFolds(group_indices, k=nfolds, list = T)
    for (fold_idx in seq_along(folds)) {
      data$fold[group_indices[folds[[fold_idx]]]] <- fold_idx
    }
  }
  data$group2 <- sub("^[^_]*_", "", data$group)
  return(data)
}
get_summ_metrics_bin <- function(y, yhat) {
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
  res <- data.frame(accuracy = accuracy,
                    F1_score = f1_score,
                    sensitivity = sensitivity,
                    specificity = specificity)
  return(res)
}
get_summ_metrics_cont <- function(y, yhat) {
  RMSE <- sqrt(mean((y-yhat)^2,na.rm=T))
  res <- data.frame(RMSE = RMSE)
  return(res)
}
cross_validate <- function(form, data, family="binomial",
                           nfolds=10, verbose=T, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  all_vars <- all.vars(formula(form))
  y_var <- all_vars[1]
  x_rand_terms <- findbars(form)
  x_group_vars <- unique(unlist(lapply(x_rand_terms, function(x) all.vars(x[[3]]))))
  x_fixed_vars <- setdiff(all_vars[-1], x_group_vars)
  groups <- c(y_var,x_group_vars)
  mix_eff <- length(groups) > 1
  for (x in x_fixed_vars) {
    if (is.factor(data[[x]])) {
      groups <- c(groups, x)
    }
  }
  data_rm <- remove_zerovar(data, groups)
  print(paste0('Removed ', nrow(data)-nrow(data_rm), ' rows with zero group-level variation in outcome'))
  data_cv <- get_folds(data_rm, groups, nfolds)
  res_kfold <- c()
  for (k in 1:nfolds) {
    if (verbose) print(paste0('Fold = ',k))
    train_df_k <- data_cv %>% filter(fold != k)
    val_df_k <- data_cv %>% filter(fold == k)
    if (mix_eff) {
      mod_k <- glmer(formula = form, data = train_df_k, family = family)
    } else {
      mod_k <- glm(formula = form, data = train_df_k, family = family)
    }
    ytrue_k <- val_df_k %>% select(y_var) %>% unlist() %>% as.integer()
    pred_k <- predict(mod_k, newdata = val_df_k, type = "response")
    if (family == "binomial") {
      yhat_k <- ifelse(pred_k>=.5,1,0) %>% as.integer()
      res_kfold <- rbind(res_kfold, get_summ_metrics_bin(ytrue_k, yhat_k))
    } else if (family == "gaussian") {
      yhat_k <- pred_k
      res_kfold <- rbind(res_kfold, get_summ_metrics_cont(ytrue_k, yhat_k))
    }
  }
  return(res_kfold)
}
predict_outcome <- function(data, nfolds = 10, family = "binomial", outcome="CQresistant", cov="Monoclonal", grp="region_grp") {
  form_fixed <- as.formula(paste0(outcome, " ~ ", cov, " + ", grp))
  form_randint <- as.formula(paste0(outcome, " ~ ", cov, " + (1|", grp,")"))
  form_randslope <- as.formula(paste0(outcome, " ~ ", cov, " + (1 + ",cov,"|", grp,")"))
  cv_res_fixed <- cross_validate(form_fixed, data, family = family, nfolds=nfolds)
  cv_res_randint <- cross_validate(form_randint, data, family = family, nfolds=nfolds)
  cv_res_randslope <- cross_validate(form_randslope, data, family = family, nfolds=nfolds)
  cv_res_fixed$Model <- "Fixed Effects"
  cv_res_randint$Model <- "Random Intercept"
  cv_res_randslope$Model <- "Random Slope"
  cv_res <- rbind(cv_res_fixed, cv_res_randint, cv_res_randslope)
  cv_res$Outcome <- outcome
  return(cv_res)
}

# Predict binary resistance outcomes with region-level
set.seed(123)
outcomes <- names(data)[endsWith(names(data), "resistant")]
cov = "Monoclonal"
grp <- "region_grp"
cv1 <- predict_outcome(data, nfolds = 25, family = "binomial", outcome = outcomes[1], cov = cov, grp = grp)
cv2 <- predict_outcome(data, nfolds = 25, family = "binomial", outcome = outcomes[2], cov = cov, grp = grp)
cv3 <- predict_outcome(data, nfolds = 25, family = "binomial", outcome = outcomes[3], cov = cov, grp = grp)
cv4 <- predict_outcome(data, nfolds = 25, family = "binomial", outcome = outcomes[4], cov = cov, grp = grp)
cv5 <- predict_outcome(data, nfolds = 25, family = "binomial", outcome = outcomes[5], cov = cov, grp = grp)
cv6 <- predict_outcome(data, nfolds = 25, family = "binomial", outcome = outcomes[6], cov = cov, grp = grp)
cv_res <- rbind(cv1, cv2, cv3, cv4, cv5, cv6)

# Plot results
#cv_res <- cv1
cv_res %>%
  ggplot(aes(y = accuracy, fill = factor(Model)))+
  facet_wrap(vars(Outcome))+
  geom_boxplot()+
  theme_bw()+
  labs(fill = "Model")+
  ylab("Accuracy")


# Misc --------------------------------------------------------------------


# Model if 3 or more resistances are observed by transformed Fws (COI) and country
set.seed(123)
nfolds <- 10
family <- "binomial"
form_fixed = CQresistant ~ Monoclonal + region_grp
form_randint = CQresistant ~ Monoclonal + (1|region_grp)
form_randslope = CQresistant ~ Monoclonal + (Monoclonal|region_grp)
cv_res_fixed <- cross_validate(form_fixed, data, family = family, nfolds=nfolds)
cv_res_randint <- cross_validate(form_randint, data, family = family, nfolds=nfolds)
cv_res_randslope <- cross_validate(form_randslope, data, family = family, nfolds=nfolds)

print(cv_res_fixed)
print(cv_res_randint)
print(cv_res_randslope)

# Make boxplots
cv_data <- rbind(cv_res_fixed, cv_res_randint, cv_res_randslope)
cv_data$Model <- rep(c("Fixed Effects Only", "Random Intercept", "Random Slope"),
                     each=nfolds)

if (family == "gaussian") {
  cv_data %>%
    ggplot(aes(y = RMSE, fill = Model))+
    geom_boxplot()+
    theme_minimal()+
    ggtitle(paste0(nfolds,'-fold CV for prediction (RMSE)'))
} else {
  cv_data %>%
    ggplot(aes(y = accuracy, fill = Model))+
    geom_boxplot()+
    theme_minimal()+
    ggtitle(paste0(nfolds,'-fold CV for prediction (accuracy)'))
  
  cv_data %>%
    ggplot(aes(y = F1_score, fill = Model))+
    geom_boxplot()+
    theme_minimal()+
    ggtitle(paste0(nfolds,'-fold CV for prediction (f1 score)'))
  

  cv_data %>%
    ggplot(aes(y = sensitivity, fill = Model))+
    geom_boxplot()+
    theme_minimal()+
    ggtitle(paste0(nfolds,'-fold CV for prediction (sensitivity)'))
  

  cv_data %>%
    ggplot(aes(y = specificity, fill = Model))+
    geom_boxplot()+
    theme_minimal()+
    ggtitle(paste0(nfolds,'-fold CV for prediction (specificity)'))
}


# Explore Effects ---------------------------------------------------------

groups <- c("CQresistant", "country")
data_rm <- remove_zerovar(data, groups)

form_fixed = CQresistant ~ Monoclonal*country
form_randint = CQresistant ~ Monoclonal + (1|country)
form_randslope = CQresistant ~ Monoclonal + (Monoclonal|country)

fullmod_fixed <- glm(formula = form_fixed, data = data, family = "binomial")
fullmod_randint <- glmer(formula = form_randint, data = data, family = "binomial")
fullmod_randslope <- glmer(formula = form_randslope, data = data, family = "binomial")


raneff <- ranef(fullmod_randslope)$country
fixed_summ <- summary(fullmod_fixed)
fixed_coef <- fixed_summ$coefficients

countries <- rownames(ranef(fullmod_randslope)$country)
randslope_ors <- exp(ranef(fullmod_randslope)$country[,2])
fixedslope_ors <- exp(as.numeric(coef(fullmod_fixed)[startsWith(names(coef(fullmod_fixed)),"Monoclonal:")]))


countries_mod <- rownames(fixed_coef)[35:65]
which(!(paste0("Monoclonal:country",countries) %in% countries_mod))


df_model_estimates <- data.frame()
df_model_estimates$countries <- rownames(raneff)
df_model_estimates$randslope <- raneff[,2]

ref_slope_ind = which(rownames(fixed_coef)=="Monoclonal")
interaction_inds <- which(startsWith(rownames(fixed_coef), "Monoclonal:country"))

countries_mod <- which(startsWith(rownames(fixed_coef), "country"))

summary(fullmod_fixed)
summary(fullmod_randint)
summary(fullmod_randslope)

VarCorr()

hist(ranef(fullmod_randslope)$country[,2])

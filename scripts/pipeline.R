library(tidyverse)
library(lme4)
library(mice)
library(caret)
library(logisticPCA)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

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
nclusters <- 5
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
  select(-c(assigned_cluster, cluster))
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

# add number of resistances
df_final <- df_clusters %>%
  mutate(nresistances = CQresistant + MQresistant + ARTresistant + PPQresistant + PYRresistant + SDXresistant) %>%
  mutate(resistant1 = ifelse(nresistances >= 1, 1, 0)) %>%
  mutate(resistant2 = ifelse(nresistances >= 2, 1, 0)) %>%
  mutate(resistant3 = ifelse(nresistances >= 3, 1, 0)) %>%
  mutate(resistant4 = ifelse(nresistances >= 4, 1, 0)) %>%
  mutate(resistant5 = ifelse(nresistances >= 5, 1, 0)) %>%
  mutate(resistant6 = ifelse(nresistances >= 6, 1, 0))

write.csv(df_final, '../data/pf7_prediction.csv')
  


# Prediction (10-fold Cross-Validation) --------------------------------------------------------------
data = df_final %>% na.omit()

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


# Model if 3 or more resistances are observed by transformed Fws (COI) and country
set.seed(123)
nfolds <- 10
family <- "binomial"
form_fixed = ARTresistant ~ Monoclonal + region_grp
form_randint = ARTresistant ~ Monoclonal + (1|region_grp)
form_randslope = ARTresistant ~ Monoclonal + (Monoclonal|region_grp)
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

groups <- c("resistant3", "country")
data_rm <- remove_zerovar(data, groups)

form_fixed = resistant3 ~ Fws_logit + country
form_randint = resistant3 ~ Fws_logit + (1|country)
form_randslope = resistant3 ~ Fws_logit + (Fws_logit|country)

fullmod_fixed <- glm(formula = form_fixed, data = data, family = "binomial")
fullmod_randint <- glmer(formula = form_randint, data = data, family = "binomial")
fullmod_randslope <- glmer(formula = form_randslope, data = data, family = "binomial")

summary(fullmod_fixed)
summary(fullmod_randint)
summary(fullmod_randslope)

hist(ranef(fullmod_randslope)$country[,2])

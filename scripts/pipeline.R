library(tidyverse)
library(lme4)

# Data Cleaning ------------------------------------------------------------
setwd('~/Dropbox/UNC/Fall2024/BIOS784/Final_Project/')
fws_data <- read.csv('pf7_fws.txt', sep = "\t")
fws_data$sample_id <- fws_data$Sample
metadata <- read.csv('Pf7-samples.csv')
df_merge <- inner_join(fws_data, metadata, by = "sample_id") %>% select(-Sample)
# df_merge <- df_merge %>%
#   mutate(across(ends_with("resistant"), ~ factor(ifelse(. == "undetermined", NA_real_, .))))
df_merge <- df_merge %>%
  mutate(across(ends_with("resistant"), ~ case_when(
    . == "sensitive" ~ 0,
    . == "resistant" ~ 1,
    . == "undetermined" ~ NA_real_
  )))
write_csv(df_merge, 'pf7_merged.csv')

# Data Processing ---------------------------------------------------------
df_model <- read.csv('pf7_merged.csv')

# logit transform 
# (deals with some of the skewness in the Fws/COI variable, which ranges from 0-1)
ep <- .001
df_model <- df_model %>%
  mutate(Fws_sc = Fws*(1-2*ep)+ep) %>%
  mutate(Fws_logit = log(Fws_sc/(1-Fws_sc))) %>%
  mutate(Monoclonal = ifelse(Fws >= .95, 1, 0))
hist(df_model$Fws)
hist((df_model$Fws_logit))

# code regions
# (can explore other region groupings)
africa <- c("Mauritania", "Gambia", "Guinea", "Kenya", "Tanzania", "Ghana", 
            "Burkina Faso", "Mali", "Malawi", "Uganda", "Democratic Republic of the Congo", 
            "Nigeria", "Madagascar", "Cameroon", "Côte d'Ivoire", "Ethiopia", 
            "Benin", "Senegal", "Gabon", "Sudan", "Mozambique")
asia <- c("Thailand", "Cambodia", "Indonesia", "Papua New Guinea", 
          "Bangladesh", "Vietnam", "Myanmar", "Laos", "India")
south_america <- c("Peru", "Colombia", "Venezuela")

df_model <- df_model %>%
  mutate(country = factor(country)) %>%
  mutate(continent = factor(case_when(country %in% africa ~ "Africa",
                                      country %in% asia ~ "Asia",
                                      TRUE ~ "South America"))) %>%
  mutate(nresistances = CQresistant + MQresistant + ARTresistant + PPQresistant + PYRresistant + SDXresistant) %>%
  mutate(resistant1 = ifelse(nresistances >= 1, 1, 0)) %>%
  mutate(resistant2 = ifelse(nresistances >= 2, 1, 0)) %>%
  mutate(resistant3 = ifelse(nresistances >= 3, 1, 0)) %>%
  mutate(resistant4 = ifelse(nresistances >= 4, 1, 0)) %>%
  mutate(resistant5 = ifelse(nresistances >= 5, 1, 0)) %>%
  mutate(resistant6 = ifelse(nresistances >= 6, 1, 0))
  
df_nona <- df_model %>% na.omit()

# Prediction (10-fold Cross-Validation) --------------------------------------------------------------
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
    probs_k <- predict(mod_k, newdata = val_df_k, type = "response")
    yhat_k <- ifelse(probs_k>=.5,1,0) %>% as.integer()
    ytrue_k <- val_df_k %>% select(y_var) %>% unlist() %>% as.integer()
    res_kfold <- rbind(res_kfold, get_summ_metrics_bin(ytrue_k, yhat_k))
  }
  return(res_kfold)
}

table(df_model$CQresistant)
table(df_model$MQresistant)
table(df_model$ARTresistant)
table(df_model$PPQresistant)
table(df_model$PYRresistant)
table(df_model$SDXresistant)
table(df_model$resistant3)

df_nona <- df_model %>% na.omit()
data = df_nona

# Model if 3 or more resistances are observed by transformed Fws (COI) and country
set.seed(123)
nfolds <- 10
form_fixed = resistant3 ~ Monoclonal + country
form_randint = resistant3 ~ Monoclonal + (1|country)
form_randslope = resistant3 ~ Monoclonal + (Monoclonal|country)
cv_res_fixed <- cross_validate(form_fixed, data, family = "binomial", nfolds=nfolds)
cv_res_randint <- cross_validate(form_randint, data, family = "binomial", nfolds=nfolds)
cv_res_randslope <- cross_validate(form_randslope, data, family = "binomial", nfolds=nfolds)

print(cv_res_fixed)
print(cv_res_randint)
print(cv_res_randslope)

# Make boxplots
cv_data <- rbind(cv_res_fixed, cv_res_randint, cv_res_randslope)
cv_data$Model <- rep(c("Fixed Effects Only", "Random Intercept", "Random Slope"),
                     each=nfolds)

cv_data %>%
  ggplot(aes(y = accuracy, fill = Model))+
  geom_boxplot()+
  theme_minimal()+
  ggtitle('10-fold CV for prediction of >=3 resistances (accuracies)')

cv_data %>%
  ggplot(aes(y = F1_score, fill = Model))+
  geom_boxplot()+
  theme_minimal()+
  ggtitle('10-fold CV for prediction of >=3 resistances (f1 scores)')


cv_data %>%
  ggplot(aes(y = sensitivity, fill = Model))+
  geom_boxplot()+
  theme_minimal()+
  ggtitle('10-fold CV for prediction of >=3 resistances (sensitivity)')


cv_data %>%
  ggplot(aes(y = specificity, fill = Model))+
  geom_boxplot()+
  theme_minimal()+
  ggtitle('10-fold CV for prediction of >=3 resistances (specificity)')



# Explore Effects ---------------------------------------------------------

groups <- c("resistant3", "country")
data_rm <- remove_zerovar(data, groups)

form_fixed = resistant3 ~ Monoclonal + country
form_randint = resistant3 ~ Monoclonal + (1|country)
form_randslope = resistant3 ~ Monoclonal + (Monoclonal|country)

fullmod_fixed <- glm(formula = form_fixed, data = data, family = "binomial")
fullmod_randint <- glmer(formula = form_randint, data = data, family = "binomial")
fullmod_randslope <- glmer(formula = form_randslope, data = data, family = "binomial")

summary(fullmod_fixed)
summary(fullmod_randint)
summary(fullmod_randslope)

hist(ranef(fullmod_randslope)$country[,2])

library(tidyverse)
library(lme4)
library(ggplot2)
library(magrittr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rlang)
library(glue)
library(patchwork)
library(purrr)

# Data Loading ------------------------------------------------------------
set.seed(123)
rm(list = ls()) # clear the environment
df_est <- read_csv("data/pf7_estimation.csv")

df_est %<>%
  mutate(region = factor(region_grp, labels = c("1", "2", "3", "4", "5")),
         country = factor(country)) %>%
  group_by(country) %>%
  filter(n() > 200)

digit_round <- 4

# Functions for estimation procedure ---------------------------------------

  ## remove countries with zero variance for outcome or exposure (doesn't make sense ow)
    ## -- I can use this function solely because the data is binary.
rem_zero_var <- function(data, res){
  data %>%
    group_by(country) %>%
    filter(!(all(Monoclonal == 0) | all(Monoclonal == 1))) %>%
    filter(!(all(!!sym(paste0(res, "resistant")) == 0 | all(!!sym(paste0(res, "resistant")) == 1)))) %>%
    ungroup()
}

estimate <- function(data, res, random, fixed_eff = TRUE){
  
  formula = case_when(fixed_eff == TRUE ~ 
                        paste(sym(paste0(res, "resistant")), 
                    "~ Monoclonal * ", sym(random)),
                      fixed_eff == FALSE ~
                        paste(sym(paste0(res, "resistant")), 
                    "~ Monoclonal + (Monoclonal|", sym(random), ")")
  )
  
  if(fixed_eff == TRUE){
    model = lm(as.formula(formula), data = data)
    estimates = coef(model)
    sigma = vcov(model)
  } else{
    model = glmer(as.formula(formula), 
                  data = data, 
                  family = binomial(link = "logit"))
    estimates = ranef(model)$country
    sigma = ranef(model, condVar = TRUE)
  }
  
  return(list(estimates, sigma, model))
}

  ## remove countries with no variance in either response or predictor
plot_est <- function(resistance, cat){
  df <- rem_zero_var(df_est, resistance)
  df$country <- droplevels(df$country)
  rand <- cat
  
    ## fixed effect
  fit <- estimate(df, resistance, rand, fixed_eff = TRUE)
  estimates <- fit[1]
  est_ref <- unname(estimates[[1]]["Monoclonal"])
  std_err_ref <- unname(sqrt(diag(summary(fit[[3]])$cov.unscaled)["Monoclonal"]))
  
  estimates <- estimates[[1]][startsWith(names(estimates[[1]]), paste0("Monoclonal:", rand))]
  estimates <- c(est_ref, unname(estimates))
  
  std_err <- sqrt(diag(summary(fit[[3]])$cov.unscaled))[startsWith(names(fit[1][[1]]), paste0("Monoclonal:", rand))]
  std_err <- c(std_err_ref, unname(std_err))
  
  df_fixed <- data.frame(a = levels(df[[rand]])[1:length(levels(df[[rand]]))], 
                        est = unname(round(estimates, digit_round)), std_err = round(std_err, digit_round))
  # df_plot %<>%
  #   mutate(CI_low = est - (1.96*std_err), 
  #          CI_up = est + (1.96*std_err)
  #   )
  
  df_fixed$method <- "non-hier"
  names(df_fixed) <- c(rand, "est", "std_err", "method") 
  
  
    ## random effects
  fit <- estimate(df, resistance, rand, fixed_eff = FALSE)
  ran_est <- fit[[1]]
  ran_sigma <- round(attr(fit[[2]]$country, "postVar")[2,2,], digit_round)
  
  ran_est <- round(ran_est$Monoclonal, digit_round)
  
  df_rand <- data.frame(a = levels(df[[rand]])[1:length(levels(df[[rand]]))], 
                        est = unname(round(ran_est, digit_round)), std_err = ran_sigma)
  
  df_rand$method <- "hier"
  names(df_rand) <- c(rand, "est", "std_err", "method") 
  
  df_plot <- rbind(df_fixed, df_rand)
  
  df_plot %<>%
    mutate(CI_low = est - (1.96*std_err), 
           CI_up = est + (1.96*std_err)
    )
  
    ## plotting the effects with the standard errors
  plot <- ggplot(df_plot, aes(x = !!sym(cat), y = est, color = method, group = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +  # Adjust position to avoid overlap
    geom_errorbar(aes(ymin = CI_low, ymax = CI_up), 
                  position = position_dodge(width = 0.5), width = 0.2) +  # Same dodge for error bars
    labs(x = cat, y = "Estimate", title = glue("Point Estimates with Confidence Intervals: {resistance} resistant")) +
    theme_minimal() +
    scale_color_manual(values = c("peachpuff4", "darkblue")) +  # Custom colors for the methods
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 15),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13))
  
  return(plot)
}


resistances <- c("CQ", "ART", "PPQ", "MQ", "PYR", "SDX")
est_plots <- list()
cat <- "country"

for(res in resistances){
  print(res)
  est_plots[[res]] = plot_est(res, cat)
  ggsave(glue("figs/{res}_est.png"), est_plots[[res]], width = 25, 
         height = 10, units = "in", dpi = 300, limitsize = FALSE)
}

#combined_plot <- reduce(est_plots, `+`) + plot_layout(ncol = 2, nrow = 1)






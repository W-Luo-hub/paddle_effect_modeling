#  HBM model fitting and comparison

library(brms)
library(dplyr)
library(loo)
library(ggplot2)

cat("Loading data...\n")
df_full <- read.csv("cleaned_data_Kirsch.csv")

df_full <- df_full %>%
  mutate(
    ExperimentName = factor(
      ExperimentName,
      levels = c("Experiment_1", "Experiment_2", "Experiment_3")
    ),
    paddle_size = factor(paddle_size, levels = c("small", "medium", "big")),
    Subject_UID = factor(Subject_UID),
    choice = as.integer(choice),
    is_medium = as.numeric(is_medium),
    is_big = as.numeric(is_big),
    is_hit = as.numeric(is_hit),
    is_exp1 = as.numeric(is_exp1),
    action_medium = as.numeric(action_medium),
    action_big = as.numeric(action_big)
  )

cat("\n--- Data check ---\n")
cat("N trials:", nrow(df_full), "\n")
cat("N subjects:", length(unique(df_full$Subject_UID)), "\n")
print(table(df_full$ExperimentName, df_full$paddle_size))
print(table(df_full$choice))
print(table(df_full$is_hit))


# MCMC settings

my_chains <- 4
my_cores  <- min(4, parallel::detectCores())
my_iter   <- 4000
my_warmup <- 2000

my_control <- list(
  adapt_delta = 0.95,
  max_treedepth = 15
)


# Priors

prior_theta_sigmap <- c(
  prior(normal(1, 0.30), class = "b", nlpar = "theta", coef = "Intercept"),
  prior(normal(0, 0.20), class = "b", nlpar = "theta"),
  prior(exponential(5), class = "sd", nlpar = "theta"),
  prior(lognormal(log(0.25), 0.50), class = "b", nlpar = "sigmap", lb = 0.01)
)

prior_with_delta <- c(
  prior_theta_sigmap,
  prior(normal(0, 0.20), class = "b", nlpar = "Delta"),
  prior(exponential(10), class = "sd", nlpar = "Delta")
)


# Model formulas

# M0: Baseline
bform_M0 <- bf(
  choice ~ Phi((x_rel - theta) / sigmap),
  theta ~ 1 + ExperimentName + (1 | Subject_UID),
  sigmap ~ 1,
  nl = TRUE,
  family = bernoulli(link = "identity")
)


# M1: Action-only
bform_M1 <- bf(
  choice ~ Phi((x_rel - (theta - Delta)) / sigmap),
  theta ~ 1 + ExperimentName + (1 | Subject_UID),
  Delta ~ 0 + action_medium + action_big +
    (0 + action_medium + action_big || Subject_UID),
  sigmap ~ 1,
  nl = TRUE,
  family = bernoulli(link = "identity")
)


# M2: Visual-hit
bform_M2 <- bf(
  choice ~ Phi((x_rel - (theta - Delta)) / sigmap),
  theta ~ 1 + ExperimentName + (1 | Subject_UID),
  Delta ~ 0 + is_hit +
    (0 + is_hit || Subject_UID),
  sigmap ~ 1,
  nl = TRUE,
  family = bernoulli(link = "identity")
)


# M3: Visual-context
bform_M3 <- bf(
  choice ~ Phi((x_rel - (theta - Delta)) / sigmap),
  theta ~ 1 + ExperimentName + (1 | Subject_UID),
  Delta ~ 0 + is_medium + is_big +
    (0 + is_medium + is_big || Subject_UID),
  sigmap ~ 1,
  nl = TRUE,
  family = bernoulli(link = "identity")
)


# M4: Full additive
bform_M4 <- bf(
  choice ~ Phi((x_rel - (theta - Delta)) / sigmap),
  theta ~ 1 + ExperimentName + (1 | Subject_UID),
  Delta ~ 0 + is_medium + is_big + is_hit + action_medium + action_big +
    (0 + is_medium + is_big + is_hit + action_medium + action_big || Subject_UID),
  sigmap ~ 1,
  nl = TRUE,
  family = bernoulli(link = "identity")
)


# Model fitting

cat("\n>>> Fitting M0: Baseline Model...\n")
fit_M0 <- brm(
  bform_M0,
  data = df_full,
  prior = prior_theta_sigmap,
  chains = my_chains,
  cores = my_cores,
  iter = my_iter,
  warmup = my_warmup,
  control = my_control,
  save_pars = save_pars(all = TRUE),
  seed = 2026,
  file = "fit_M0",
  file_refit = "always"
)

cat("\n>>> Fitting M1: Action-only Model...\n")
fit_M1 <- brm(
  bform_M1,
  data = df_full,
  prior = prior_with_delta,
  chains = my_chains,
  cores = my_cores,
  iter = my_iter,
  warmup = my_warmup,
  control = my_control,
  save_pars = save_pars(all = TRUE),
  seed = 2026,
  file = "fit_M1",
  file_refit = "always"
)

cat("\n>>> Fitting M2: Visual-hit Model...\n")
fit_M2 <- brm(
  bform_M2,
  data = df_full,
  prior = prior_with_delta,
  chains = my_chains,
  cores = my_cores,
  iter = my_iter,
  warmup = my_warmup,
  control = my_control,
  save_pars = save_pars(all = TRUE),
  seed = 2026,
  file = "fit_M2",
  file_refit = "always"
)

cat("\n>>> Fitting M3: Visual-context Model...\n")
fit_M3 <- brm(
  bform_M3,
  data = df_full,
  prior = prior_with_delta,
  chains = my_chains,
  cores = my_cores,
  iter = my_iter,
  warmup = my_warmup,
  control = my_control,
  save_pars = save_pars(all = TRUE),
  seed = 2026,
  file = "fit_M3",
  file_refit = "always"
)

cat("\n>>> Fitting M4: Full Additive Model...\n")
fit_M4 <- brm(
  bform_M4,
  data = df_full,
  prior = prior_with_delta,
  chains = my_chains,
  cores = my_cores,
  iter = my_iter,
  warmup = my_warmup,
  control = my_control,
  save_pars = save_pars(all = TRUE),
  seed = 2026,
  file = "fit_M4",
  file_refit = "always"
)


# Save models

cat("\nSaving all models...\n")
save(
  fit_M0, fit_M1, fit_M2, fit_M3, fit_M4,
  file = "Fitted_Models_HBM_Psychometric.RData"
)
cat("Saved to 'Fitted_Models_HBM_Psychometric.RData'.\n")


# Posterior summaries


cat("\n>>> [M0] Baseline Model:\n")
print(summary(fit_M0))

cat("\n>>> [M1] Action-only Model:\n")
print(summary(fit_M1))

cat("\n>>> [M2] Visual-hit Model:\n")
print(summary(fit_M2))

cat("\n>>> [M3] Visual-context Model:\n")
print(summary(fit_M3))

cat("\n>>> [M4] Full Additive Model:\n")
print(summary(fit_M4))


# LOO and WAIC


cat("\n>>> Computing PSIS-LOO...\n")
fit_M0 <- add_criterion(fit_M0, "loo")
fit_M1 <- add_criterion(fit_M1, "loo")
fit_M2 <- add_criterion(fit_M2, "loo")
fit_M3 <- add_criterion(fit_M3, "loo")
fit_M4 <- add_criterion(fit_M4, "loo")

cat("\n>>> Computing WAIC...\n")
fit_M0 <- add_criterion(fit_M0, "waic")
fit_M1 <- add_criterion(fit_M1, "waic")
fit_M2 <- add_criterion(fit_M2, "waic")
fit_M3 <- add_criterion(fit_M3, "waic")
fit_M4 <- add_criterion(fit_M4, "waic")


# Pareto k diagnostics

cat("\n--- Pareto k diagnostics ---\n")

cat("\nM0:\n")
print(pareto_k_table(fit_M0$criteria$loo))

cat("\nM1:\n")
print(pareto_k_table(fit_M1$criteria$loo))

cat("\nM2:\n")
print(pareto_k_table(fit_M2$criteria$loo))

cat("\nM3:\n")
print(pareto_k_table(fit_M3$criteria$loo))

cat("\nM4:\n")
print(pareto_k_table(fit_M4$criteria$loo))


# Model comparison

cat("\n--- LOO comparison ---\n")
loo_results <- loo_compare(
  fit_M0, fit_M1, fit_M2, fit_M3, fit_M4,
  criterion = "loo"
)
print(loo_results)

cat("\n--- WAIC comparison ---\n")
waic_results <- loo_compare(
  fit_M0, fit_M1, fit_M2, fit_M3, fit_M4,
  criterion = "waic"
)
print(waic_results)

write.csv(as.data.frame(loo_results), "LOO_model_comparison.csv")
write.csv(as.data.frame(waic_results), "WAIC_model_comparison.csv")

cat("\nLOO and WAIC comparison saved as CSV.\n")


# ELPD interpretation

summarize_elpd <- function(comparison_results, metric_name = "LOO") {
  
  cat(sprintf("\n[ %s ELPD difference summary ]\n", metric_name))
  
  res_df <- as.data.frame(comparison_results)
  res_df$model <- rownames(res_df)
  best_model <- res_df$model[1]
  
  cat(sprintf("Best predictive model: %s\n\n", best_model))
  
  for (i in 2:nrow(res_df)) {
    
    comp_model <- res_df$model[i]
    diff_val <- abs(res_df$elpd_diff[i])
    se_val <- res_df$se_diff[i]
    ratio <- diff_val / se_val
    
    cat(sprintf("-> %s vs %s:\n", best_model, comp_model))
    cat(sprintf(
      "   |Delta ELPD| = %.2f, SE = %.2f, ratio = %.2f\n",
      diff_val, se_val, ratio
    ))
    
    if (ratio >= 2) {
      cat("   Interpretation: relatively stable predictive advantage.\n")
    } else if (ratio >= 1) {
      cat("   Interpretation: possible but uncertain predictive advantage.\n")
    } else {
      cat("   Interpretation: weak or inconclusive predictive difference.\n")
    }
  }
}

summarize_elpd(loo_results, "PSIS-LOO")
summarize_elpd(waic_results, "WAIC")


# Final save

save(
  fit_M0, fit_M1, fit_M2, fit_M3, fit_M4,
  loo_results, waic_results,
  file = "Fitted_Models_HBM_Psychometric_WithCriteria.RData"
)



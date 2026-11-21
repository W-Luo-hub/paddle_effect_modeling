# Figure 4

library(brms)
library(dplyr)
library(ggplot2)
library(posterior)

proj_root <- if (file.exists("cleaned_data_Kirsch.csv")) {
  "."
} else if (file.exists("../cleaned_data_Kirsch.csv")) {
  ".."
} else {
  "."
}

cat("Loading models and data...\n")

rdata <- file.path(proj_root, "Fitted_Models_HBM_Psychometric_WithCriteria.RData")
fit_rds <- file.path(proj_root, "fit_M4.rds")

if (file.exists(rdata)) {
  load(rdata)
  if (!exists("fit_M4")) {
    stop("Loaded ", rdata, " but object fit_M4 was not found.", sep = "")
  }
} else if (file.exists(fit_rds)) {
  fit_M4 <- readRDS(fit_rds)
} else {
  stop(
    "Neither Fitted_Models_HBM_Psychometric_WithCriteria.RData nor fit_M4.rds ",
    "was found. Run Main_Analysis_HBM.R first."
  )
}

df_full <- read.csv(file.path(proj_root, "cleaned_data_Kirsch.csv")) %>%
  mutate(
    ExperimentName = factor(
      ExperimentName,
      levels = c("Experiment_1", "Experiment_2", "Experiment_3")
    ),
    paddle_size = factor(
      paddle_size,
      levels = c("small", "medium", "big")
    ),
    Subject_UID = factor(Subject_UID),
    choice = as.integer(choice),
    is_medium = as.numeric(is_medium),
    is_big = as.numeric(is_big),
    is_hit = as.numeric(is_hit),
    is_exp1 = as.numeric(is_exp1),
    action_medium = as.numeric(action_medium),
    action_big = as.numeric(action_big)
  )

theme_set(
  theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
)

pal_paddle <- c(Small = "#3D6399", Medium = "#DD8452", Big = "#D94738")

save_tiff <- function(filename, plot, width, height, dpi = 600) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "tiff",
    compression = "lzw"
  )
}

cat("Plotting Figure 4: subject x paddle predicted vs observed...\n")

make_subject_paddle_ppc_m4 <- function(fit, data, ndraws = 500) {

  ep <- posterior_epred(
    fit,
    newdata = data,
    ndraws = ndraws,
    re_formula = NULL
  )

  stopifnot(ncol(ep) == nrow(data))

  data_idx <- data %>%
    mutate(.row_id = row_number())

  group_info <- data_idx %>%
    group_by(ExperimentName, Subject_UID, paddle_size) %>%
    summarise(
      obs = mean(choice),
      n = n(),
      rows = list(.row_id),
      .groups = "drop"
    )

  pred_mean <- vapply(group_info$rows, function(idx) {
    draw_means <- rowMeans(ep[, idx, drop = FALSE])
    mean(draw_means)
  }, numeric(1))

  bind_cols(
    group_info %>% select(-rows),
    pred = pred_mean
  )
}

ppc_M4_subject <- make_subject_paddle_ppc_m4(fit_M4, df_full, ndraws = 500) %>%
  mutate(
    paddle_size = recode(
      as.character(paddle_size),
      "small" = "Small",
      "medium" = "Medium",
      "big" = "Big"
    ),
    paddle_size = factor(paddle_size, levels = c("Small", "Medium", "Big")),
    ExperimentName = recode(
      as.character(ExperimentName),
      "Experiment_1" = "Experiment 1",
      "Experiment_2" = "Experiment 2",
      "Experiment_3" = "Experiment 3"
    ),
    ExperimentName = factor(
      ExperimentName,
      levels = c("Experiment 1", "Experiment 2", "Experiment 3")
    )
  )

cat("\n--- Figure 4 data checks ---\n")
cat("Total points:", nrow(ppc_M4_subject), "(expected 108)\n")
n_by_exp <- ppc_M4_subject %>% count(ExperimentName)
print(n_by_exp)
cat("Trials per point:\n")
print(summary(ppc_M4_subject$n))
cat("Any NA in obs/pred:", anyNA(ppc_M4_subject$obs) || anyNA(ppc_M4_subject$pred), "\n")
cat(
  "Observed range: [",
  round(min(ppc_M4_subject$obs), 3), ", ",
  round(max(ppc_M4_subject$obs), 3), "]\n",
  sep = ""
)
cat(
  "Predicted range: [",
  round(min(ppc_M4_subject$pred), 3), ", ",
  round(max(ppc_M4_subject$pred), 3), "]\n",
  sep = ""
)

stopifnot(nrow(ppc_M4_subject) == 108)
stopifnot(all(n_by_exp$n == 36))
stopifnot(all(ppc_M4_subject$n == 135))
stopifnot(!anyNA(ppc_M4_subject$obs), !anyNA(ppc_M4_subject$pred))
stopifnot(all(ppc_M4_subject$pred >= 0 & ppc_M4_subject$pred <= 1))
stopifnot(all(ppc_M4_subject$obs >= 0 & ppc_M4_subject$obs <= 1))

fit_stats <- function(d) {
  data.frame(
    n = nrow(d),
    r = cor(d$obs, d$pred, use = "complete.obs"),
    RMSE = sqrt(mean((d$pred - d$obs)^2))
  )
}

cat("\n--- Fit statistics (console only) ---\n")
cat("Overall:\n")
print(fit_stats(ppc_M4_subject))
cat("By experiment:\n")
print(
  ppc_M4_subject %>%
    group_by(ExperimentName) %>%
    group_modify(~ fit_stats(.x)) %>%
    ungroup()
)

fig4 <- ggplot(
  ppc_M4_subject,
  aes(x = obs, y = pred, color = paddle_size)
) +
  annotate(
    "segment",
    x = 0.35, y = 0.35,
    xend = 0.85, yend = 0.85,
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.85
  ) +
  facet_wrap(~ ExperimentName, nrow = 1) +
  scale_x_continuous(
    limits = c(0.35, 0.85),
    breaks = seq(0.40, 0.80, by = 0.10),
    labels = function(x) sprintf("%.2f", x)
  ) +
  scale_y_continuous(
    limits = c(0.35, 0.85),
    breaks = seq(0.40, 0.80, by = 0.10),
    labels = function(x) sprintf("%.2f", x)
  ) +
  scale_color_manual(values = pal_paddle) +
  labs(
    x = "Observed P(test judged faster)",
    y = "Predicted P(test judged faster)",
    color = "Paddle size"
  ) +
  theme(
    aspect.ratio = 1,
    panel.spacing.x = unit(1.2, "lines"),
    plot.margin = margin(t = 2, r = 5.5, b = 2, l = 5.5, unit = "pt"),
    legend.position = c(0.99, 0.06),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.85), color = NA),
    legend.key = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.85), color = NA)
  )

save_tiff(
  filename = file.path(proj_root, "Figure4_M4_Predicted_vs_Observed.tiff"),
  plot = fig4,
  width = 10,
  height = 5
)

cat(
  "\nSaved: ",
  file.path(proj_root, "Figure4_M4_Predicted_vs_Observed.tiff"),
  "\n",
  sep = ""
)

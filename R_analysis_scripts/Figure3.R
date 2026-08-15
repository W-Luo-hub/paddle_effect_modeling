# Figure 3

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)

proj_root <- if (file.exists("cleaned_data_Kirsch.csv")) {
  "."
} else if (file.exists("../cleaned_data_Kirsch.csv")) {
  ".."
} else {
  "."
}

out_dir <- file.path(proj_root, "figure3_standalone")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading models and data...\n")

rdata <- file.path(proj_root, "Fitted_Models_HBM_Psychometric_WithCriteria.RData")
if (!file.exists(rdata)) {
  stop("File not found: ", rdata, ". Run 03_Main_Analysis_HBM.R first.", sep = "")
}
load(rdata)

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

cat("Plotting Figure 3: observed vs predicted psychometric curves...\n")

make_group_ppc_m4 <- function(fit, data, ndraws = 500) {

  ep <- posterior_epred(
    fit,
    newdata = data,
    ndraws = ndraws,
    re_formula = NULL
  )

  data_idx <- data %>%
    mutate(.row_id = row_number())

  group_info <- data_idx %>%
    group_by(ExperimentName, paddle_size, x_rel, v_test) %>%
    summarise(
      obs = mean(choice),
      n = n(),
      rows = list(.row_id),
      .groups = "drop"
    )

  pred_stats <- lapply(group_info$rows, function(idx) {
    draw_means <- rowMeans(ep[, idx, drop = FALSE])
    c(
      pred_mean = mean(draw_means),
      pred_lower = unname(quantile(draw_means, 0.025)),
      pred_upper = unname(quantile(draw_means, 0.975))
    )
  })

  pred_stats <- as.data.frame(do.call(rbind, pred_stats))

  bind_cols(
    group_info %>% select(-rows),
    pred_stats
  )
}

ppc_M4_main <- make_group_ppc_m4(fit_M4, df_full, ndraws = 500) %>%
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
    )
  )

fig3 <- ggplot(
  ppc_M4_main,
  aes(x = v_test, color = paddle_size, fill = paddle_size)
) +
  geom_ribbon(
    aes(ymin = pred_lower, ymax = pred_upper),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(
    aes(y = pred_mean),
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = obs),
    size = 1.8,
    alpha = 0.85
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 35,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  facet_wrap(~ ExperimentName, nrow = 1) +
  scale_x_continuous(
    breaks = c(7, 14, 21, 28, 35, 42, 49, 56, 63)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  scale_color_manual(values = pal_paddle) +
  scale_fill_manual(values = pal_paddle) +
  labs(
    x = "Test speed (cm/s)",
    y = "P(test judged faster)",
    color = "Paddle size",
    fill = "Paddle size"
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
  filename = file.path(out_dir, "Figure3_M4_Psychometric_Curves.tiff"),
  plot = fig3,
  width = 10,
  height = 5
)



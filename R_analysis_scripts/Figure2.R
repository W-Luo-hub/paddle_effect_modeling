# Figure 2

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(ggdist)

proj_root <- if (file.exists("cleaned_data_Kirsch.csv")) {
  "."
} else if (file.exists("../cleaned_data_Kirsch.csv")) {
  ".."
} else {
  "."
}

out_dir <- file.path(proj_root, "figure2_standalone")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading models and data...\n")

rdata <- file.path(proj_root, "Fitted_Models_HBM_Psychometric_WithCriteria.RData")
if (!file.exists(rdata)) {
  stop("File not found: ", rdata, ". Run 03_Main_Analysis_HBM.R first.", sep = "")
}
load(rdata)

theme_set(
  theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      panel.grid.major.y = element_blank()
    )
)

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

# Parameter map
param_map <- tibble(
  b_col = c(
    "b_Delta_is_medium",
    "b_Delta_is_big",
    "b_Delta_is_hit",
    "b_Delta_action_medium",
    "b_Delta_action_big"
  ),
  parameter = c(
    "Visual-context (M)",
    "Visual-context (B)",
    "Visual-hit",
    "Action-only (M)",
    "Action-only (B)"
  ),
  re_term = c(
    "is_medium",
    "is_big",
    "is_hit",
    "action_medium",
    "action_big"
  ),
  re_idx = 1:5
)

param_levels <- rev(param_map$parameter)

col_slab_fill <- "#3D6399"
col_slab_line <- "#2C3E50"
col_subj <- "#2C3E50"

pick_re_list <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  if ("Subject_UID" %in% names(x)) {
    out <- x$Subject_UID
    if (is.list(out) && "Delta" %in% names(out)) return(out$Delta)
    return(out)
  }
  if ("Delta" %in% names(x)) return(x$Delta)
  x
}

match_re_term_name <- function(nm, param_map) {
  if (is.null(nm) || is.na(nm) || nm == "") return(NA_character_)
  for (j in seq_len(nrow(param_map))) {
    rt <- param_map$re_term[j]
    if (nm %in% c(rt, paste0("Delta_", rt), param_map$b_col[j], as.character(j))) {
      return(rt)
    }
    if (grepl(rt, nm, fixed = TRUE)) return(rt)
    if (grepl(paste0("Delta_", rt), nm, fixed = TRUE)) return(rt)
  }
  idx <- suppressWarnings(as.integer(nm))
  if (!is.na(idx) && idx >= 1L && idx <= nrow(param_map)) {
    return(param_map$re_term[idx])
  }
  NA_character_
}

#' Normalize brms RE structures to list(term -> matrix/data.frame)
normalize_re_list <- function(re_obj, param_map) {
  if (is.null(re_obj) || !length(re_obj)) return(list())

  if (is.array(re_obj) && length(dim(re_obj)) == 3L) {
    dn1 <- dimnames(re_obj)[[1]]
    dn2 <- dimnames(re_obj)[[2]]
    dn3 <- dimnames(re_obj)[[3]]
    est_idx <- if (!is.null(dn2)) which(dn2 == "Estimate")[1] else 1L
    if (is.na(est_idx)) est_idx <- 1L
    out <- list()
    for (k in seq_len(dim(re_obj)[3L])) {
      key <- if (!is.null(dn3) && length(dn3) >= k) {
        match_re_term_name(dn3[k], param_map)
      } else {
        param_map$re_term[k]
      }
      if (is.na(key)) next
      vals <- re_obj[, est_idx, k]
      out[[key]] <- data.frame(
        Estimate = as.numeric(vals),
        row.names = if (!is.null(dn1)) dn1 else NULL
      )
    }
    return(out)
  }

  if (is.matrix(re_obj) && !is.null(colnames(re_obj)) && ncol(re_obj) == nrow(param_map)) {
    out <- list()
    for (k in seq_len(ncol(re_obj))) {
      key <- match_re_term_name(colnames(re_obj)[k], param_map)
      if (is.na(key)) key <- param_map$re_term[k]
      out[[key]] <- data.frame(
        Estimate = re_obj[, k],
        row.names = rownames(re_obj)
      )
    }
    return(out)
  }

  if (is.data.frame(re_obj)) {
    if ("Estimate" %in% names(re_obj) && any(c("level", "Level", "group", "Group") %in% names(re_obj))) {
      lvl_col <- intersect(c("level", "Level"), names(re_obj))[1]
      par_col <- intersect(c("term", "Term", "coef", "Coef", "parameter", "Parameter"), names(re_obj))[1]
      if (!is.na(lvl_col) && !is.na(par_col)) {
        out <- split(re_obj, re_obj[[par_col]])
        names(out) <- vapply(names(out), match_re_term_name, character(1), param_map = param_map)
        out <- out[!is.na(names(out))]
        names(out) <- names(out)
        return(out)
      }
    }
    if ("Estimate" %in% names(re_obj) && !is.null(rownames(re_obj))) {
      return(list(`__single_df__` = re_obj))
    }
  }

  if (!is.list(re_obj)) {
    return(list())
  }

  out <- list()
  nms <- names(re_obj)
  for (i in seq_along(re_obj)) {
    nm <- if (length(nms) >= i) nms[i] else as.character(i)
    key <- match_re_term_name(nm, param_map)
    if (is.na(key) && length(re_obj) == nrow(param_map)) {
      key <- param_map$re_term[i]
    }
    if (!is.na(key)) out[[key]] <- re_obj[[i]]
  }
  out
}

rows_from_term_table <- function(el, parameter, levels_vec) {
  if (is.data.frame(el) && "Estimate" %in% names(el)) {
    ids <- rownames(el)
    if (is.null(ids) || !length(ids)) return(NULL)
    return(tibble(
      Subject_UID = ids,
      subj_mean = as.numeric(el$Estimate),
      parameter = factor(parameter, levels = levels_vec)
    ))
  }
  NULL
}

rows_from_term_matrix <- function(mat, pop, parameter, levels_vec) {
  mat <- as.matrix(mat)
  n_draws <- length(pop)
  if (nrow(mat) == n_draws) {
    ids <- colnames(mat)
    sm <- colMeans(pop + mat)
  } else if (ncol(mat) == n_draws) {
    ids <- rownames(mat)
    sm <- rowMeans(pop + mat)
  } else {
    return(NULL)
  }
  if (is.null(ids) || !length(ids)) return(NULL)
  tibble(
    Subject_UID = ids,
    subj_mean = as.numeric(sm),
    parameter = factor(parameter, levels = levels_vec)
  )
}

#' Population-level posterior draws (long format)
get_population_draws <- function(fit, param_map, levels_vec) {
  d <- as.data.frame(as_draws_df(fit))
  missing <- setdiff(param_map$b_col, names(d))
  if (length(missing)) {
    stop("Missing posterior columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  d %>%
    select(all_of(param_map$b_col)) %>%
    pivot_longer(everything(), names_to = "b_col", values_to = "value") %>%
    left_join(param_map %>% select(b_col, parameter), by = "b_col") %>%
    mutate(parameter = factor(parameter, levels = levels_vec))
}

#' Subject-level coefficients from coef/ranef
get_subject_means <- function(fit, param_map, levels_vec) {
  rows <- list()

  # method 1: coef(summary=TRUE)
  cf <- tryCatch(
    coef(fit, nlpar = "Delta", group = "Subject_UID", summary = TRUE),
    error = function(e) NULL
  )
  cf_obj <- pick_re_list(cf)
  if (!is.null(cf_obj) && length(cf_obj)) {
    cat("coef(summary=TRUE) raw names: ", paste(head(names(cf_obj), 8), collapse = ", "), "\n", sep = "")
    cf_norm <- normalize_re_list(cf_obj, param_map)
    cat("coef aligned: ", paste(names(cf_norm), collapse = ", "), "\n", sep = "")
    for (i in seq_len(nrow(param_map))) {
      term <- param_map$re_term[i]
      el <- cf_norm[[term]]
      if (is.null(el)) next
      r <- rows_from_term_table(el, param_map$parameter[i], levels_vec)
      if (!is.null(r)) rows[[length(rows) + 1]] <- r
    }
  }

  if (length(rows)) {
    return(bind_rows(rows))
  }

  # method 1b: ranef(summary=TRUE) + fixef
  fe <- tryCatch(fixef(fit, nlpar = "Delta", summary = TRUE), error = function(e) NULL)
  rf <- tryCatch(ranef(fit, nlpar = "Delta", summary = TRUE), error = function(e) NULL)
  rf_obj <- pick_re_list(rf)
  if (!is.null(fe) && !is.null(rf_obj) && length(rf_obj)) {
    rf_norm <- normalize_re_list(rf_obj, param_map)
    for (i in seq_len(nrow(param_map))) {
      term <- param_map$re_term[i]
      el <- rf_norm[[term]]
      if (is.null(el) || !"Estimate" %in% names(el)) next
      bname <- param_map$b_col[i]
      pop_est <- if (bname %in% rownames(fe)) fe[bname, "Estimate"] else NA_real_
      if (is.na(pop_est)) {
        hit <- grep(term, rownames(fe), value = TRUE)[1]
        if (!is.na(hit)) pop_est <- fe[hit, "Estimate"]
      }
      if (is.na(pop_est)) next
      rows[[length(rows) + 1]] <- tibble(
        Subject_UID = rownames(el),
        subj_mean = pop_est + el$Estimate,
        parameter = factor(param_map$parameter[i], levels = levels_vec)
      )
    }
  }

  if (length(rows)) {
    return(bind_rows(rows))
  }

  # method 2: ranef(FALSE) + population fixed effects
  draws <- as.data.frame(as_draws_df(fit))
  re <- tryCatch(ranef(fit, nlpar = "Delta", summary = FALSE), error = function(e) NULL)
  re_obj <- if (!is.null(re)) pick_re_list(re) else NULL
  if (!is.null(re_obj) && length(re_obj)) {
    cat("ranef(FALSE) raw names: ", paste(head(names(re_obj), 8), collapse = ", "), "\n", sep = "")
    re_norm <- normalize_re_list(re_obj, param_map)
    cat("ranef aligned: ", paste(names(re_norm), collapse = ", "), "\n", sep = "")
    for (i in seq_len(nrow(param_map))) {
      term <- param_map$re_term[i]
      el <- re_norm[[term]]
      if (is.null(el)) next
      if (is.data.frame(el) && "Estimate" %in% names(el)) {
        r <- rows_from_term_table(el, param_map$parameter[i], levels_vec)
      } else {
        r <- rows_from_term_matrix(el, draws[[param_map$b_col[i]]], param_map$parameter[i], levels_vec)
      }
      if (!is.null(r)) rows[[length(rows) + 1]] <- r
    }
  }

  if (length(rows)) {
    return(bind_rows(rows))
  }

  # method 3: posterior columns r_Subject_UID*
  all_re <- grep("^r_Subject_UID", names(draws), value = TRUE)
  delta_re <- all_re[!grepl("__theta", all_re, fixed = TRUE)]
  cat("Posterior Delta RE columns: ", length(delta_re), "\n", sep = "")
  if (length(delta_re)) {
    for (col in delta_re) {
      for (i in seq_len(nrow(param_map))) {
        if (!grepl(param_map$re_term[i], col, fixed = TRUE) &&
            !grepl(paste0("Delta_", param_map$re_term[i]), col, fixed = TRUE)) {
          next
        }
        inner <- sub(".*\\[|\\]$", "", col)
        parts <- strsplit(inner, ",", fixed = TRUE)[[1]]
        subj <- parts[grepl("^Experiment", parts)][1]
        if (is.na(subj)) subj <- parts[length(parts)]
        rows[[length(rows) + 1]] <- tibble(
          Subject_UID = subj,
          subj_mean = mean(draws[[param_map$b_col[i]]] + draws[[col]]),
          parameter = factor(param_map$parameter[i], levels = levels_vec)
        )
      }
    }
    if (length(rows)) {
      return(bind_rows(rows) %>% distinct(Subject_UID, parameter, .keep_all = TRUE))
    }
  }

  tibble()
}

# Figure 2

cat("Plotting Figure 2: M4 Delta raincloud...\n")

draws_long <- get_population_draws(fit_M4, param_map, param_levels)
subject_means <- get_subject_means(fit_M4, param_map, param_levels)
show_subjects <- nrow(subject_means) > 0

if (!show_subjects) {
  warning("Could not extract subject-level Delta; plotting population posterior only.", call. = FALSE)
} else {
  cat(
    "Subject-level posterior means: ", nrow(subject_means), " rows (",
    n_distinct(subject_means$Subject_UID), " subjects)\n",
    sep = ""
  )
  set.seed(2026)
  subject_means <- subject_means %>%
    mutate(
      y_base = as.numeric(parameter) - 0.22,
      y_dot = y_base + runif(n(), min = -0.06, max = 0.06)
    )
}

pop_summary <- draws_long %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value),
    lo = quantile(value, 0.025),
    hi = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    label_ci = sprintf("%.2f [%.2f, %.2f]", mean, lo, hi)
  )

x_hi <- max(c(draws_long$value, pop_summary$hi), na.rm = TRUE)
x_lo <- min(c(draws_long$value, pop_summary$lo), na.rm = TRUE)
x_pad <- max(0.12, diff(range(c(x_lo, x_hi))) * 0.28)
x_right <- x_hi + x_pad

y_expand_bottom <- 0.42
y_expand_top <- 0.48
n_param <- length(param_levels)
y_axis_lo <- 1 - y_expand_bottom
y_axis_hi <- n_param + y_expand_top

zero_line <- geom_segment(
  x = 0,
  xend = 0,
  y = y_axis_lo,
  yend = y_axis_hi,
  inherit.aes = FALSE,
  linetype = "dashed",
  linewidth = 0.4,
  color = "grey50"
)

fig2 <- ggplot(draws_long, aes(x = value, y = parameter)) +
  stat_halfeye(
    fill = col_slab_fill,
    color = col_slab_line,
    slab_fill = col_slab_fill,
    slab_alpha = 0.32,
    slab_color = NA,
    point_color = col_slab_line,
    point_size = 1.35,
    interval_color = col_slab_line,
    interval_size = 0.5,
    .width = 0.95,
    side = "top",
    height = 0.4
  )

if (show_subjects) {
  fig2 <- fig2 +
    geom_point(
      data = subject_means,
      aes(x = subj_mean, y = y_dot),
      inherit.aes = FALSE,
      size = 0.95,
      alpha = 0.45,
      color = col_subj,
      stroke = 0.2
    )
}

fig2 <- fig2 +
  zero_line +
  geom_text(
    data = pop_summary,
    aes(x = Inf, y = parameter, label = label_ci),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0.5,
    size = 2.85,
    color = "grey20"
  ) +
  annotate(
    "text",
    x = Inf,
    y = n_param + 0.35,
    label = "Estimate [95% CrI]",
    hjust = 1,
    vjust = 0,
    size = 3.1,
    fontface = "bold",
    color = "grey25"
  ) +
  scale_y_discrete(expand = expansion(add = c(y_expand_bottom, y_expand_top))) +
  coord_cartesian(
    xlim = c(x_lo - 0.02, x_right),
    clip = "off"
  ) +
  labs(
    x = "Latent PSE shift",
    y = NULL
  ) +
  theme(
    plot.margin = margin(t = 28, r = 16, b = 28, l = 5.5, unit = "pt"),
    axis.text.y = element_text(size = 10, color = "black")
  )


save_tiff(
  filename = file.path(out_dir, "Figure2_M4_Posterior_Intervals.tiff"),
  plot = fig2,
  width = 8.5,
  height = 5
)



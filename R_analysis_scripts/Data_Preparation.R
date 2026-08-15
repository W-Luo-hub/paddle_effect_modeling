#  Data preparation

library(dplyr)
library(haven)

cat("Reading raw data...\n")
df_raw_full <- read_sav("raw_data_Kirsch&Kunde.sav")

df_full <- df_raw_full %>%
  mutate(
    ExperimentName = as.character(zap_labels(ExperimentName)),
    paddle_size = as.character(zap_labels(paddle_size)),
    
    x_rel = as.numeric(test_V_rel_to_stadard_V),
    
    v_std = 35,
    v_test = 35 * x_rel,
    
    choice = as.integer(as.numeric(zap_labels(judgment))),
    
    is_hit = as.integer(as.numeric(zap_labels(hit))),
    
    is_exp1 = ifelse(ExperimentName == "Experiment_1", 1, 0),
    
    paddle_size = factor(paddle_size, levels = c("small", "medium", "big")),
    
    ExperimentName = factor(
      ExperimentName,
      levels = c("Experiment_1", "Experiment_2", "Experiment_3")
    ),
    
    Subject_UID = factor(paste(ExperimentName, Subject, sep = "_")),
    
    is_medium = ifelse(paddle_size == "medium", 1, 0),
    is_big = ifelse(paddle_size == "big", 1, 0),
    
    action_medium = is_medium * is_exp1,
    action_big = is_big * is_exp1
  ) %>%
  filter(
    !is.na(choice),
    !is.na(x_rel),
    !is.na(paddle_size),
    !is.na(is_hit)
  ) %>%
  select(
    ExperimentName, Subject, Subject_UID, Block, Trial,
    paddle_size, is_medium, is_big,
    hit, is_hit, is_exp1, action_medium, action_big,
    x_rel, v_std, v_test, choice
  )

cat("\n--- Data structure check ---\n")
cat("N subjects:", length(unique(df_full$Subject_UID)), "\n")
cat("N trials:", nrow(df_full), "\n\n")

cat("--- Key variable check ---\n")
print(table(df_full$ExperimentName, df_full$paddle_size))
print(table(df_full$choice))
print(table(df_full$is_hit))

cat("\n--- First 6 rows ---\n")
print(head(df_full))

cat("\n--- Column types ---\n")
str(df_full)

write.csv(df_full, "cleaned_data_Kirsch.csv", row.names = FALSE)



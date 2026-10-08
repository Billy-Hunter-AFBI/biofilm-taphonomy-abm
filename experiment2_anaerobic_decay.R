
# ------------------------------------------------------------
# Time to anoxia
#
# NA indicates that the simulation did not reach the
# provisional anoxia threshold during the model run.
# ------------------------------------------------------------

summary_time <- aggregate(

  results["time_to_anoxia"],

  by = results[group_columns],

  FUN = function(x) {

    if (all(is.na(x))) {
      return(NA_real_)
    }

    mean(x, na.rm = TRUE)
  },

  na.action = na.pass
)

names(summary_time)[
  names(summary_time) == "time_to_anoxia"
] <- "mean_time_to_anoxia"


# ------------------------------------------------------------
# Proportion of simulations reaching anoxia
# ------------------------------------------------------------

summary_anoxia_frequency <- aggregate(

  results["time_to_anoxia"],

  by = results[group_columns],

  FUN = function(x) {
    mean(!is.na(x))
  },

  na.action = na.pass
)

names(summary_anoxia_frequency)[
  names(summary_anoxia_frequency) == "time_to_anoxia"
] <- "proportion_reaching_anoxia"


# ============================================================
# 10. Combine summary tables
# ============================================================

summary_results <- Reduce(

  function(x, y) {

    merge(
      x,
      y,
      by = group_columns,
      all = TRUE
    )
  },

  list(
    summary_mean,
    summary_sd,
    summary_se,
    summary_time,
    summary_anoxia_frequency
  )
)


# ============================================================
# 11. Save parameter-space summary
# ============================================================

write.csv(

  summary_results,

  file.path(
    out_dir,
    "experiment_2_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 12. Sensitivity relative to complete anaerobic arrest
# ============================================================
#
# Compare each anaerobic decay treatment with r = 0.
#
# Delta survival:
#
#   P(r) - P(r = 0)
#
# Negative values indicate reduced preservation relative
# to complete anaerobic arrest.
#
# ============================================================

baseline <- subset(
  summary_results,
  anaerobic_decay_ratio == 0,
  select = c(
    "D_O2",
    "k_resp_aero",
    "final_survival"
  )
)

names(baseline)[3] <-
  "baseline_survival"


sensitivity_results <- merge(

  summary_results,

  baseline,

  by = c(
    "D_O2",
    "k_resp_aero"
  ),

  all.x = TRUE
)


sensitivity_results$delta_survival <-
  sensitivity_results$final_survival -
  sensitivity_results$baseline_survival


# ============================================================
# 13. Save sensitivity results
# ============================================================

write.csv(

  sensitivity_results,

  file.path(
    out_dir,
    "experiment_2_sensitivity.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 14. Overall preservation sensitivity
# ============================================================
#
# Average preservation across the full environmental
# parameter grid for each anaerobic decay ratio.
#
# This is a descriptive summary, not an estimate of
# preservation probability in natural environments.
# ============================================================

overall_sensitivity <- aggregate(

  final_survival ~ anaerobic_decay_ratio,

  data = summary_results,

  FUN = mean
)

names(overall_sensitivity)[2] <-
  "mean_survival"


# ------------------------------------------------------------
# Range across environmental parameter combinations
# ------------------------------------------------------------

survival_min <- aggregate(

  final_survival ~ anaerobic_decay_ratio,

  data = summary_results,

  FUN = min
)

names(survival_min)[2] <-
  "minimum_survival"


survival_max <- aggregate(

  final_survival ~ anaerobic_decay_ratio,

  data = summary_results,

  FUN = max
)

names(survival_max)[2] <-
  "maximum_survival"


overall_sensitivity <- Reduce(

  function(x, y) {

    merge(
      x,
      y,
      by = "anaerobic_decay_ratio"
    )
  },

  list(
    overall_sensitivity,
    survival_min,
    survival_max
  )
)


# ============================================================
# 15. Save overall sensitivity
# ============================================================

write.csv(

  overall_sensitivity,

  file.path(
    out_dir,
    "experiment_2_overall_sensitivity.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 16. Basic diagnostics
# ============================================================

cat(
  "\n\nEXPERIMENT 2 SUMMARY\n",
  "====================\n\n"
)

cat(
  "Mean survival by anaerobic decay ratio:\n"
)

print(overall_sensitivity)


cat(
  "\nOverall survival range:\n"
)

print(
  range(
    results$final_survival,
    na.rm = TRUE
  )
)


cat(
  "\nSensitivity relative to complete arrest:\n"
)

print(
  aggregate(
    delta_survival ~ anaerobic_decay_ratio,
    data = sensitivity_results,
    FUN = mean
  )
)


# ============================================================
# 17. Reproducibility information
# ============================================================

capture.output(

  sessionInfo(),

  file = file.path(
    out_dir,
    "experiment_2_session_info.txt"
  )
)


# Save experimental design

write.csv(

  parameter_grid,

  file.path(
    out_dir,
    "experiment_2_parameter_grid.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 18. Completion
# ============================================================

cat(
  "\nExperiment 2 complete.\n",
  "Results saved to:\n",
  out_dir,
  "\n",
  sep = ""
)

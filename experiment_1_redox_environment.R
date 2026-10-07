# ============================================================
# EXPERIMENT 1
# REDOX PHASE SPACE AND TISSUE PRESERVATION
# ============================================================
#
# Question:
#
# How does preservation emerge from the balance between
# oxygen resupply and aerobic oxygen consumption?
#
# Experimental factors:
#
#   1. Oxygen diffusivity / resupply (D_O2)
#   2. Aerobic respiration rate (k_resp_aero)
#
# Primary response:
#
#   Final whole-organism tissue survival
#
# Secondary diagnostics:
#
#   - cumulative anoxia
#   - final anoxia
#   - time to anoxia
#   - final biofilm
#   - tissue-specific survival
#
# ============================================================


# ============================================================
# 1. Setup
# ============================================================

base_dir <- getwd()

fig_dir <- file.path(
  base_dir,
  "figures",
  "experiment_1"
)

out_dir <- file.path(
  base_dir,
  "outputs",
  "experiment_1"
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


library(ggplot2)


# ------------------------------------------------------------
# Load body plan and shared model
# ------------------------------------------------------------

source("hallucigenia_bodyplan.R")
source("taphonomy_model_core.R")


# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

set.seed(42)


# ============================================================
# 2. Construct model organism
# ============================================================

bp <- make_hallucigenia_offset(
  N = 100,
  target_row = 10
)


# ============================================================
# 3. Experimental design
# ============================================================
#
# IMPORTANT:
#
# These values represent dimensionless model parameter space.
# They should not currently be interpreted as empirically
# calibrated oxygen diffusion or respiration rates.
#
# The ranges are deliberately broad enough to test whether
# different redox regimes emerge.
#
# ============================================================

D_O2_values <- seq(
  0.05,
  0.50,
  length.out = 7
)

resp_values <- seq(
  0.025,
  0.100,
  length.out = 7
)

n_replicates <- 20


parameter_grid <- expand.grid(

  D_O2 = D_O2_values,

  k_resp_aero = resp_values,

  replicate = seq_len(n_replicates),

  KEEP.OUT.ATTRS = FALSE
)


cat(
  "\nExperiment 1\n",
  "------------\n",
  "Oxygen-resupply levels: ",
  length(D_O2_values),
  "\n",
  "Respiration levels: ",
  length(resp_values),
  "\n",
  "Replicates: ",
  n_replicates,
  "\n",
  "Total simulations: ",
  nrow(parameter_grid),
  "\n\n",
  sep = ""
)


# ============================================================
# 4. Storage
# ============================================================

results_list <- vector(
  "list",
  nrow(parameter_grid)
)


# ============================================================
# 5. Run simulations
# ============================================================

for (i in seq_len(nrow(parameter_grid))) {

  pars <- parameter_grid[i, ]


  # ----------------------------------------------------------
  # Progress indicator
  # ----------------------------------------------------------

  if (
    i == 1 ||
    i %% 50 == 0 ||
    i == nrow(parameter_grid)
  ) {

    cat(
      "Running simulation",
      i,
      "of",
      nrow(parameter_grid),
      "\n"
    )
  }


  # ----------------------------------------------------------
  # Run model
  # ----------------------------------------------------------

  sim <- run_taphonomy_model(

    bp = bp,

    steps = 120,

    D_C = 0.01,

    D_O2 = pars$D_O2,

    bc_top = 1,

    k_resp_aero = pars$k_resp_aero,

    k_resp_anaer = 0.005,

    k_sed = 0.002,

    o2_crit = 0.2,

    anoxia_sharpness = 20,

    k_bio = 0.03,

    biofilm_protection = 0.7,

    anaerobic_decay_ratio = 0.25,

    phi_min = 0.2,

    phi_max = 1,

    record_history = TRUE
  )


  # ----------------------------------------------------------
  # Tissue-specific results
  # ----------------------------------------------------------

  tissue_results <- as.list(
    sim$final_tissue_survival
  )


  # ----------------------------------------------------------
  # Store simulation
  # ----------------------------------------------------------

  results_list[[i]] <- c(

    list(

      D_O2 = pars$D_O2,

      k_resp_aero = pars$k_resp_aero,

      replicate = pars$replicate,

      final_survival =
        sim$final_survival,

      cumulative_anoxia =
        sim$cumulative_anoxia,

      final_anoxia =
        sim$final_anoxia,

      final_o2 =
        sim$final_o2,

      time_to_anoxia =
        sim$time_to_anoxia,

      final_biofilm =
        sim$final_biofilm
    ),

    tissue_results
  )
}


# ============================================================
# 6. Combine results
# ============================================================

results <- do.call(
  rbind.data.frame,
  results_list
)

rownames(results) <- NULL


# Ensure numeric columns remain numeric after binding

numeric_columns <- setdiff(
  names(results),
  character(0)
)

results[numeric_columns] <- lapply(
  results[numeric_columns],
  as.numeric
)


# ============================================================
# 7. Save raw simulation results
# ============================================================

write.csv(
  results,
  file.path(
    out_dir,
    "experiment_1_raw_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. Summarise parameter combinations
# ============================================================

summary_results <- aggregate(

  cbind(
    final_survival,
    cumulative_anoxia,
    final_anoxia,
    final_o2,
    final_biofilm
  ) ~

    D_O2 +
    k_resp_aero,

  data = results,

  FUN = mean
)


# ------------------------------------------------------------
# Standard deviation of preservation
# ------------------------------------------------------------

survival_sd <- aggregate(

  final_survival ~

    D_O2 +
    k_resp_aero,

  data = results,

  FUN = sd
)

names(survival_sd)[3] <-
  "final_survival_sd"


# ------------------------------------------------------------
# Standard error
# ------------------------------------------------------------

survival_se <- aggregate(

  final_survival ~

    D_O2 +
    k_resp_aero,

  data = results,

  FUN = function(x) {

    sd(x) / sqrt(length(x))
  }
)

names(survival_se)[3] <-
  "final_survival_se"


# ------------------------------------------------------------
# Mean time to anoxia
#
# NA values mean the simulation never crossed the provisional
# mean-anoxia threshold of 0.5.
# ------------------------------------------------------------

time_summary <- aggregate(

  time_to_anoxia ~

    D_O2 +
    k_resp_aero,

  data = results,

  FUN = function(x) {

    if (all(is.na(x))) {

      NA_real_

    } else {

      mean(
        x,
        na.rm = TRUE
      )
    }
  },

  na.action = na.pass
)


names(time_summary)[3] <-
  "mean_time_to_anoxia"


# ------------------------------------------------------------
# Probability of reaching predominantly anoxic conditions
#
# This is NOT preservation probability.
#
# It is simply the proportion of replicates in which mean
# anoxia over the original body area reached >= 0.5.
# ------------------------------------------------------------

anoxia_probability <- aggregate(

  time_to_anoxia ~

    D_O2 +
    k_resp_aero,

  data = results,

  FUN = function(x) {

    mean(
      !is.na(x)
    )
  },

  na.action = na.pass
)


names(anoxia_probability)[3] <-
  "proportion_reaching_anoxia"


# ------------------------------------------------------------
# Merge summaries
# ------------------------------------------------------------

summary_results <- merge(
  summary_results,
  survival_sd,
  by = c(
    "D_O2",
    "k_resp_aero"
  )
)

summary_results <- merge(
  summary_results,
  survival_se,
  by = c(
    "D_O2",
    "k_resp_aero"
  )
)

summary_results <- merge(
  summary_results,
  time_summary,
  by = c(
    "D_O2",
    "k_resp_aero"
  )
)

summary_results <- merge(
  summary_results,
  anoxia_probability,
  by = c(
    "D_O2",
    "k_resp_aero"
  )
)


# ============================================================
# 9. Save summary table
# ============================================================

write.csv(
  summary_results,
  file.path(
    out_dir,
    "experiment_1_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. Figure theme
# ============================================================

theme_experiment <- theme_bw(
  base_size = 12
) +

  theme(

    panel.grid.minor =
      element_blank(),

    panel.grid.major =
      element_blank(),

    axis.title =
      element_text(
        size = 12
      ),

    axis.text =
      element_text(
        size = 10
      ),

    plot.title =
      element_text(
        size = 14,
        face = "bold"
      ),

    plot.subtitle =
      element_text(
        size = 11
      ),

    legend.title =
      element_text(
        size = 11
      )
  )


# ============================================================
# 11. FIGURE 1
# Preservation phase space
# ============================================================
#
# This is the principal figure for Experiment 1.
#
# If the hypothesis is correct, we expect preservation to
# increase toward conditions combining:
#
#   low oxygen resupply
#   high oxygen demand
#
# ============================================================

p_preservation <- ggplot(

  summary_results,

  aes(
    x = D_O2,
    y = k_resp_aero,
    fill = final_survival
  )

) +

  geom_tile() +

  geom_contour(

    aes(
      z = final_survival
    ),

    colour = "black",

    bins = 8,

    linewidth = 0.35
  ) +

  scale_fill_viridis_c(

    name =
      "Final tissue\nsurvival",

    limits = c(
      0,
      1
    )
  ) +

  labs(

    title =
      "Preservation across redox parameter space",

    subtitle =
      "Mean final tissue survival across replicate simulations",

    x =
      "Oxygen resupply (D[O2])",

    y =
      "Aerobic oxygen consumption (k[resp])"
  ) +

  theme_experiment


ggsave(

  filename = file.path(
    fig_dir,
    "figure_1_preservation_phase_space.png"
  ),

  plot = p_preservation,

  width = 7,

  height = 6,

  dpi = 600
)


ggsave(

  filename = file.path(
    fig_dir,
    "figure_1_preservation_phase_space.pdf"
  ),

  plot = p_preservation,

  width = 7,

  height = 6
)


# ============================================================
# 12. FIGURE 2
# Cumulative anoxia phase space
# ============================================================

p_anoxia <- ggplot(

  summary_results,

  aes(
    x = D_O2,
    y = k_resp_aero,
    fill = cumulative_anoxia
  )

) +

  geom_tile() +

  geom_contour(

    aes(
      z = cumulative_anoxia
    ),

    colour = "black",

    bins = 8,

    linewidth = 0.35
  ) +

  scale_fill_viridis_c(

    name =
      "Cumulative\nanoxia",

    limits = c(
      0,
      1
    )
  ) +

  labs(

    title =
      "Emergent anoxia across redox parameter space",

    subtitle =
      "Mean anoxic exposure over the original body area",

    x =
      "Oxygen resupply (D[O2])",

    y =
      "Aerobic oxygen consumption (k[resp])"
  ) +

  theme_experiment


ggsave(

  filename = file.path(
    fig_dir,
    "figure_2_cumulative_anoxia.png"
  ),

  plot = p_anoxia,

  width = 7,

  height = 6,

  dpi = 600
)


ggsave(

  filename = file.path(
    fig_dir,
    "figure_2_cumulative_anoxia.pdf"
  ),

  plot = p_anoxia,

  width = 7,

  height = 6
)


# ============================================================
# 13. FIGURE 3
# Preservation versus cumulative anoxia
# ============================================================
#
# This figure asks whether the emergent preservation response
# is actually associated with the redox history generated by
# the model.
#
# Individual replicate simulations are shown.
# ============================================================

p_mechanism <- ggplot(

  results,

  aes(
    x = cumulative_anoxia,
    y = final_survival
  )

) +

  geom_point(
    alpha = 0.25,
    size = 1.5
  ) +

  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 0.8
  ) +

  labs(

    title =
      "Relationship between anoxic exposure and preservation",

    subtitle =
      "Each point represents one simulation",

    x =
      "Cumulative anoxic exposure",

    y =
      "Final tissue survival"
  ) +

  coord_cartesian(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +

  theme_experiment


ggsave(

  filename = file.path(
    fig_dir,
    "figure_3_anoxia_vs_preservation.png"
  ),

  plot = p_mechanism,

  width = 7,

  height = 6,

  dpi = 600
)


ggsave(

  filename = file.path(
    fig_dir,
    "figure_3_anoxia_vs_preservation.pdf"
  ),

  plot = p_mechanism,

  width = 7,

  height = 6
)


# ============================================================
# 14. FIGURE 4
# Stochastic variability in preservation
# ============================================================
#
# Useful diagnostic:
#
# Does the transition zone also show increased variability
# among replicate simulations?
#
# ============================================================

p_variability <- ggplot(

  summary_results,

  aes(
    x = D_O2,
    y = k_resp_aero,
    fill = final_survival_sd
  )

) +

  geom_tile() +

  scale_fill_viridis_c(

    name =
      "SD of tissue\nsurvival"
  ) +

  labs(

    title =
      "Variability in preservation outcomes",

    subtitle =
      "Standard deviation among replicate simulations",

    x =
      "Oxygen resupply (D[O2])",

    y =
      "Aerobic oxygen consumption (k[resp])"
  ) +

  theme_experiment


ggsave(

  filename = file.path(
    fig_dir,
    "figure_4_preservation_variability.png"
  ),

  plot = p_variability,

  width = 7,

  height = 6,

  dpi = 600
)


# ============================================================
# 15. Tissue-specific summary
# ============================================================

tissue_columns <- grep(
  "^tissue_",
  names(results),
  value = TRUE
)


tissue_summary_list <- lapply(

  tissue_columns,

  function(tissue_name) {

    tmp <- aggregate(

      results[[tissue_name]],

      by = list(

        D_O2 =
          results$D_O2,

        k_resp_aero =
          results$k_resp_aero
      ),

      FUN = mean
    )

    names(tmp)[3] <-
      "survival"

    tmp$tissue <-
      tissue_name

    tmp
  }
)


tissue_summary <- do.call(
  rbind,
  tissue_summary_list
)


write.csv(

  tissue_summary,

  file.path(
    out_dir,
    "experiment_1_tissue_specific_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 16. FIGURE 5
# Tissue-specific preservation
# ============================================================
#
# This is initially a diagnostic / supplementary figure.
#
# It tells us whether environmental regime alters not simply
# how much tissue survives, but which tissue classes survive.
# ============================================================

p_tissues <- ggplot(

  tissue_summary,

  aes(
    x = D_O2,
    y = k_resp_aero,
    fill = survival
  )

) +

  geom_tile() +

  facet_wrap(
    ~ tissue
  ) +

  scale_fill_viridis_c(

    name =
      "Tissue\nsurvival",

    limits = c(
      0,
      1
    )
  ) +

  labs(

    title =
      "Tissue-specific preservation across redox space",

    x =
      "Oxygen resupply (D[O2])",

    y =
      "Aerobic oxygen consumption (k[resp])"
  ) +

  theme_experiment


ggsave(

  filename = file.path(
    fig_dir,
    "figure_5_tissue_specific_preservation.png"
  ),

  plot = p_tissues,

  width = 9,

  height = 7,

  dpi = 600
)


ggsave(

  filename = file.path(
    fig_dir,
    "figure_5_tissue_specific_preservation.pdf"
  ),

  plot = p_tissues,

  width = 9,

  height = 7
)


# ============================================================
# 17. Basic diagnostic statistics
# ============================================================

cat(
  "\n\nEXPERIMENT 1 SUMMARY\n",
  "====================\n\n"
)

cat(
  "Final survival range:\n"
)

print(
  range(
    summary_results$final_survival,
    na.rm = TRUE
  )
)


cat(
  "\nCumulative anoxia range:\n"
)

print(
  range(
    summary_results$cumulative_anoxia,
    na.rm = TRUE
  )
)


cat(
  "\nCorrelation between cumulative anoxia ",
  "and final survival:\n",
  sep = ""
)

print(
  cor(
    results$cumulative_anoxia,
    results$final_survival,
    use = "complete.obs"
  )
)


cat(
  "\nMaximum between-replicate SD:\n"
)

print(
  max(
    summary_results$final_survival_sd,
    na.rm = TRUE
  )
)


# ============================================================
# 18. Save session information
# ============================================================

capture.output(

  sessionInfo(),

  file = file.path(
    out_dir,
    "experiment_1_session_info.txt"
  )
)


# ============================================================
# End
# ============================================================

cat(
  "\nExperiment 1 complete.\n",
  "Results written to: ",
  out_dir,
  "\n",
  "Figures written to: ",
  fig_dir,
  "\n",
  sep = ""
)
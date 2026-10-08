
# ============================================================
# FIGURE 2
# ANAEROBIC DECAY SENSITIVITY
# SIX-PANEL PRESERVATION LANDSCAPE
# ============================================================
#
# Input:
#   outputs/experiment_2/experiment_2_summary.csv
#
# Outputs:
#   figures/experiment_2/
#     figure_2_anaerobic_decay_sensitivity.png
#     figure_2_anaerobic_decay_sensitivity.pdf
#
# ============================================================


# ============================================================
# 1. SETUP
# ============================================================

library(ggplot2)

base_dir <- getwd()

input_file <- file.path(
  base_dir,
  "outputs",
  "experiment_2",
  "experiment_2_summary.csv"
)

fig_dir <- file.path(
  base_dir,
  "figures",
  "experiment_2"
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Experiment 2 summary file not found: ",
    input_file
  )
}


# ============================================================
# 2. LOAD DATA
# ============================================================

results <- read.csv(
  input_file
)


# ============================================================
# 3. VALIDATE DATA
# ============================================================

required_columns <- c(
  "D_O2",
  "k_resp_aero",
  "anaerobic_decay_ratio",
  "final_survival"
)

missing_columns <- setdiff(
  required_columns,
  names(results)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyNA(results[required_columns])) {
  stop("Required plotting variables contain missing values.")
}

if (any(
  results$final_survival < 0 |
  results$final_survival > 1
)) {
  stop("Final survival values must lie between 0 and 1.")
}


# ============================================================
# 4. PREPARE FACET LABELS
# ============================================================

decay_levels <- c(
  0, 0.10, 0.25, 0.50, 0.75, 1.00
)

if (!setequal(
  unique(results$anaerobic_decay_ratio),
  decay_levels
)) {
  stop("Unexpected or missing anaerobic decay treatments.")
}

results$decay_panel <- factor(

  results$anaerobic_decay_ratio,

  levels = decay_levels,

  labels = c(
    "r = 0.00",
    "r = 0.10",
    "r = 0.25",
    "r = 0.50",
    "r = 0.75",
    "r = 1.00"
  )
)


# ============================================================
# 5. VALIDATE PARAMETER GRID
# ============================================================

expected_D <- sort(unique(results$D_O2))
expected_R <- sort(unique(results$k_resp_aero))

expected_grid <- expand.grid(
  D_O2 = expected_D,
  k_resp_aero = expected_R,
  anaerobic_decay_ratio = decay_levels
)

actual_keys <- paste(
  results$D_O2,
  results$k_resp_aero,
  results$anaerobic_decay_ratio
)

expected_keys <- paste(
  expected_grid$D_O2,
  expected_grid$k_resp_aero,
  expected_grid$anaerobic_decay_ratio
)

if (
  nrow(results) != nrow(expected_grid) ||
  anyDuplicated(actual_keys) ||
  !setequal(actual_keys, expected_keys)
) {
  stop("Incomplete or duplicated parameter grid.")
}


# ============================================================
# 6. FIGURE 2
# ============================================================
#
# Each panel shows mean final tissue survival.
#
# X-axis:
#   Oxygen resupply
#
# Y-axis:
#   Aerobic oxygen consumption
#
# Colour:
#   Final tissue survival, fixed from 0 to 1.
#
# ============================================================

figure_2 <- ggplot(

  results,

  aes(
    x = D_O2,
    y = k_resp_aero,
    fill = final_survival
  )

) +

  # ----------------------------------------------------------
  # Preservation landscape
  # ----------------------------------------------------------

  geom_tile() +

  # ----------------------------------------------------------
  # Contours of equal preservation
  # ----------------------------------------------------------

  geom_contour(

    aes(z = final_survival),

    breaks = seq(
      0.1,
      0.9,
      by = 0.1
    ),

    colour = "black",

    linewidth = 0.3,

    alpha = 0.65
  ) +

  # ----------------------------------------------------------
  # Six anaerobic decay treatments
  # ----------------------------------------------------------

  facet_wrap(
    ~ decay_panel,
    ncol = 3
  ) +

  # ----------------------------------------------------------
  # Common colour scale
  # ----------------------------------------------------------

  scale_fill_viridis_c(

    name = "Final tissue\nsurvival",

    limits = c(0, 1),

    breaks = seq(
      0,
      1,
      by = 0.2
    ),

    option = "D"
  ) +

  # ----------------------------------------------------------
  # Axis labels
  # ----------------------------------------------------------

  labs(

    x = expression(
      "Oxygen resupply (" * D[O[2]] * ")"
    ),

    y = expression(
      "Aerobic oxygen consumption (" * k[resp] * ")"
    )
  ) +

  # ----------------------------------------------------------
  # Theme
  # ----------------------------------------------------------

  theme_classic(
    base_size = 12
  ) +

  theme(

    strip.background = element_blank(),

    strip.text = element_text(
      size = 12,
      face = "bold"
    ),

    panel.spacing = grid::unit(
      1.1,
      "lines"
    ),

    axis.text = element_text(
      size = 10,
      colour = "black"
    ),

    axis.title = element_text(
      size = 12
    ),

    legend.title = element_text(
      size = 11
    ),

    legend.text = element_text(
      size = 10
    ),

    legend.key.height = grid::unit(
      2.5,
      "cm"
    ),

    plot.margin = margin(
      10, 15, 10, 10
    )
  )


# ============================================================
# 7. DISPLAY FIGURE
# ============================================================

print(figure_2)


# ============================================================
# 8. SAVE PNG
# ============================================================

ggsave(

  filename = file.path(
    fig_dir,
    "figure_2_anaerobic_decay_sensitivity.png"
  ),

  plot = figure_2,

  width = 10,

  height = 7,

  units = "in",

  dpi = 600,

  bg = "white"
)


# ============================================================
# 9. SAVE VECTOR PDF
# ============================================================

ggsave(

  filename = file.path(
    fig_dir,
    "figure_2_anaerobic_decay_sensitivity.pdf"
  ),

  plot = figure_2,

  width = 10,

  height = 7,

  units = "in",

  bg = "white"
)


# ============================================================
# 10. DIAGNOSTIC SUMMARY
# ============================================================

cat(
  "\nFIGURE 2: ANAEROBIC DECAY SENSITIVITY\n",
  "====================================\n\n"
)

cat("Final survival range by treatment:\n")

print(

  aggregate(

    final_survival ~ anaerobic_decay_ratio,

    data = results,

    FUN = function(x) {
      c(
        min = min(x),
        max = max(x),
        mean = mean(x)
      )
    }
  )
)


cat(
  "\nFigure 2 saved to:\n",
  fig_dir,
  "\n",
  sep = ""
)

# ============================================================
# END
# ============================================================

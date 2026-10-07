# ============================================================
# FIGURE 1
# PRESERVATION ACROSS REDOX PARAMETER SPACE
# ============================================================
#
# Experiment 1:
# Interaction between oxygen resupply and aerobic oxygen
# consumption in determining final tissue preservation.
#
# Input:
#   outputs/experiment_1/experiment_1_summary.csv
#
# Output:
#   figures/experiment_1/
#       figure_1_preservation_phase_space.png
#       figure_1_preservation_phase_space.pdf
#
# ============================================================


# ============================================================
# 1. Setup
# ============================================================

library(ggplot2)

base_dir <- getwd()

input_file <- file.path(
  base_dir,
  "outputs",
  "experiment_1",
  "experiment_1_summary.csv"
)

fig_dir <- file.path(
  base_dir,
  "figures",
  "experiment_1"
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Load Experiment 1 summary
# ============================================================

results <- read.csv(
  input_file
)


# ============================================================
# 3. Basic checks
# ============================================================

required_columns <- c(
  "D_O2",
  "k_resp_aero",
  "final_survival"
)

missing_columns <- setdiff(
  required_columns,
  names(results)
)

if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing required columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 4. Inspect preservation range
# ============================================================

cat(
  "\nFinal tissue survival range:\n"
)

print(
  range(
    results$final_survival,
    na.rm = TRUE
  )
)


# ============================================================
# 5. Figure 1
# ============================================================
#
# Tile colour:
#   Mean final tissue survival
#
# Contours:
#   Lines of equal tissue survival
#
# Preservation is kept on an absolute 0-1 scale so that
# subsequent experiments can be directly compared.
#
# ============================================================

figure_1 <- ggplot(

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
  # Preservation contours
  # ----------------------------------------------------------

  geom_contour(

    aes(
      z = final_survival
    ),

    colour = "black",

    bins = 8,

    linewidth = 0.35,

    alpha = 0.7
  ) +


  # ----------------------------------------------------------
  # Colour scale
  # ----------------------------------------------------------

  scale_fill_viridis_c(

    name =
      "Final tissue\nsurvival",

    limits = c(
      0,
      1
    ),

    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +


  # ----------------------------------------------------------
  # Axis labels
  # ----------------------------------------------------------

  labs(

    x =
      expression(
        "Oxygen resupply (" *
        D[O[2]] *
        ")"
      ),

    y =
      expression(
        "Aerobic oxygen consumption (" *
        k[resp] *
        ")"
      )
  ) +


  # ----------------------------------------------------------
  # Theme
  # ----------------------------------------------------------

  theme_classic(
    base_size = 12
  ) +

  theme(

    axis.title.x =
      element_text(
        margin = margin(
          t = 10
        )
      ),

    axis.title.y =
      element_text(
        margin = margin(
          r = 10
        )
      ),

    axis.text =
      element_text(
        colour = "black"
      ),

    legend.title =
      element_text(
        size = 11
      ),

    legend.text =
      element_text(
        size = 10
      ),

    legend.key.height =
      grid::unit(
        1.5,
        "cm"
      ),

    plot.margin =
      margin(
        10,
        15,
        10,
        10
      )
  )


# ============================================================
# 6. Display figure
# ============================================================

print(
  figure_1
)


# ============================================================
# 7. Save publication-resolution PNG
# ============================================================

ggsave(

  filename = file.path(
    fig_dir,
    "figure_1_preservation_phase_space.png"
  ),

  plot = figure_1,

  width = 7,

  height = 5.5,

  units = "in",

  dpi = 600
)


# ============================================================
# 8. Save vector PDF
# ============================================================

ggsave(

  filename = file.path(
    fig_dir,
    "figure_1_preservation_phase_space.pdf"
  ),

  plot = figure_1,

  width = 7,

  height = 5.5,

  units = "in"
)


# ============================================================
# End
# ============================================================

cat(
  "\nFigure 1 saved to:\n",
  fig_dir,
  "\n"
)
# ============================================================
# TAPHONOMY MODEL CORE
# ============================================================
#
# Shared mechanistic functions for Burgess Shale-style
# taphonomic simulations.
#
# Experimental parameter grids, plotting, statistical analysis
# and scenario definitions should NOT be placed in this file.
#
# Requires:
#   hallucigenia_bodyplan.R
#
# ============================================================


# ============================================================
# 1. Spatial transport functions
# ============================================================

# ------------------------------------------------------------
# Neumann (zero-flux) Laplacian
#
# Used for transport of tissue-derived organic carbon.
# Material cannot leave the model domain through its edges.
# ------------------------------------------------------------

laplacian_neumann_c <- function(Z) {

  N <- nrow(Z)
  M <- ncol(Z)

  Z_up    <- rbind(Z[2, ], Z[-N, ])
  Z_down  <- rbind(Z[-1, ], Z[N - 1, ])
  Z_right <- cbind(Z[, -1], Z[, M - 1])
  Z_left  <- cbind(Z[, 2], Z[, -M])

  Z_up + Z_down + Z_right + Z_left - 4 * Z
}


# ------------------------------------------------------------
# Mixed-boundary Laplacian for oxygen
#
# Top boundary:
#   fixed oxygen concentration (Dirichlet boundary)
#
# Other boundaries:
#   zero flux (Neumann boundary)
#
# This represents oxygen resupply from the overlying water.
# ------------------------------------------------------------

laplacian_mixed_o2 <- function(Z, bc_top = 1) {

  N <- nrow(Z)
  M <- ncol(Z)

  Z_up    <- rbind(rep(bc_top, M), Z[-N, ])
  Z_down  <- rbind(Z[-1, ], Z[N - 1, ])
  Z_right <- cbind(Z[, -1], Z[, M - 1])
  Z_left  <- cbind(Z[, 2], Z[, -M])

  Z_up + Z_down + Z_right + Z_left - 4 * Z
}


# ------------------------------------------------------------
# Oxygen diffusion
#
# Sub-stepping is used when D_O2 becomes large to improve
# numerical stability.
# ------------------------------------------------------------

diffuse_o2 <- function(O2,
                       D_O2,
                       bc_top = 1) {

  n_sub <- max(1L, ceiling(D_O2 / 0.20))
  D_sub <- D_O2 / n_sub

  for (s in seq_len(n_sub)) {

    O2 <- O2 +
      D_sub * laplacian_mixed_o2(O2, bc_top)

    O2 <- pmax(pmin(O2, 1), 0)

    # Reinstate fixed oxygen concentration at sediment surface
    O2[1, ] <- bc_top
  }

  O2
}


# ============================================================
# 2. Redox functions
# ============================================================

# ------------------------------------------------------------
# Continuous anoxia response
#
# Converts oxygen concentration into an anoxia weighting
# between 0 and 1.
#
# This avoids imposing a hard oxic/anoxic threshold.
# ------------------------------------------------------------

anoxia_field <- function(O2,
                         o2_crit = 0.2,
                         sharpness = 20) {

  1 / (1 + exp(sharpness * (O2 - o2_crit)))
}


# ============================================================
# 3. Body-plan preparation
# ============================================================

# ------------------------------------------------------------
# Position Hallucigenia within the model domain
#
# The existing Hallucigenia body plan is shifted toward the
# sediment-water interface.
# ------------------------------------------------------------

make_hallucigenia_offset <- function(N = 100,
                                     target_row = 10) {

  bp <- make_hallucigenia(N)

  # Original body plan is centred approximately on row 50
  shift <- target_row - 50

  shift_matrix <- function(mat, fill = 0) {

    result <- matrix(fill, N, ncol(mat))

    for (r in seq_len(N)) {

      r_new <- r + shift

      if (r_new >= 1 && r_new <= N) {
        result[r_new, ] <- mat[r, ]
      }
    }

    result
  }

  bp$C <- shift_matrix(
    bp$C,
    fill = 0
  )

  bp$k_tissue <- shift_matrix(
    bp$k_tissue,
    fill = 0
  )

  bp$tissue <- shift_matrix(
    bp$tissue,
    fill = 0
  )

  bp$tissue_labels <- shift_matrix(
    bp$tissue_labels,
    fill = "empty"
  )

  bp
}


# ============================================================
# 4. Preservation metrics
# ============================================================

# ------------------------------------------------------------
# Whole-organism tissue survival
#
# IMPORTANT:
# Preservation is deliberately calculated independently of
# oxygen/anoxia.
#
# This avoids defining successful preservation partly by the
# environmental mechanism hypothesised to cause preservation.
# ------------------------------------------------------------

compute_tissue_survival <- function(C,
                                    C_initial,
                                    body_mask) {

  sum(C[body_mask]) /
    (sum(C_initial[body_mask]) + 1e-12)
}


# ------------------------------------------------------------
# Tissue-specific survival
#
# Returns survival for each original tissue class.
# ------------------------------------------------------------

compute_tissue_specific_survival <- function(C,
                                             C_initial,
                                             tissue_ref) {

  tissue_ids <- sort(
    unique(tissue_ref[tissue_ref > 0])
  )

  survival <- sapply(
    tissue_ids,
    function(k) {

      mask <- tissue_ref == k

      sum(C[mask]) /
        (sum(C_initial[mask]) + 1e-12)
    }
  )

  names(survival) <- paste0(
    "tissue_", tissue_ids
  )

  survival
}


# ============================================================
# 5. Main simulation
# ============================================================

run_taphonomy_model <- function(

    bp,

    # Simulation duration
    steps = 120,

    # Organic-carbon transport
    D_C = 0.01,

    # Oxygen transport
    D_O2 = 0.10,
    bc_top = 1,

    # Respiration
    k_resp_aero = 0.015,
    k_resp_anaer = 0.005,
    k_sed = 0.002,

    # Redox response
    o2_crit = 0.2,
    anoxia_sharpness = 20,

    # Biofilm
    k_bio = 0.03,
    biofilm_protection = 0.7,

    # Relative tissue degradation under anoxia
    #
    # 0 = complete arrest
    # 1 = same intrinsic degradation rate as oxic tissue
    anaerobic_decay_ratio = 0.25,

    # Random spatial heterogeneity
    phi_min = 0.2,
    phi_max = 1,

    # Diagnostics
    record_history = TRUE

) {

  # ----------------------------------------------------------
  # Initial state
  # ----------------------------------------------------------

  C_initial <- bp$C
  C         <- bp$C

  k_tissue  <- bp$k_tissue
  tissue_ref <- bp$tissue

  body_mask <- C_initial > 0

  N <- nrow(C)

  # Oxygenated starting environment
  O2 <- matrix(
    1,
    nrow = N,
    ncol = N
  )

  # Spatial heterogeneity in microbial/biofilm potential
  phi <- matrix(
    runif(
      N * N,
      min = phi_min,
      max = phi_max
    ),
    nrow = N,
    ncol = N
  )

  # No initial biofilm
  B <- matrix(
    0,
    nrow = N,
    ncol = N
  )

  # Fixed oxygenated upper boundary
  O2[1, ] <- bc_top


  # ----------------------------------------------------------
  # Diagnostic storage
  # ----------------------------------------------------------

  mean_anoxia <- numeric(steps)
  mean_o2     <- numeric(steps)
  mean_biofilm <- numeric(steps)
  total_survival <- numeric(steps)

  tissue_ids <- sort(
    unique(tissue_ref[tissue_ref > 0])
  )

  tissue_survival <- matrix(
    NA_real_,
    nrow = steps,
    ncol = length(tissue_ids)
  )

  colnames(tissue_survival) <- paste0(
    "tissue_", tissue_ids
  )


  # ==========================================================
  # Simulation loop
  # ==========================================================

  for (t in seq_len(steps)) {

    # --------------------------------------------------------
    # Current redox state
    # --------------------------------------------------------

    anox <- anoxia_field(
      O2,
      o2_crit = o2_crit,
      sharpness = anoxia_sharpness
    )


    # --------------------------------------------------------
    # Transport of tissue-derived organic material
    # --------------------------------------------------------

    C <- C +
      D_C * laplacian_neumann_c(C)

    C <- pmax(
      pmin(C, 1),
      0
    )


    # --------------------------------------------------------
    # Oxygen resupply
    # --------------------------------------------------------

    O2 <- diffuse_o2(
      O2,
      D_O2 = D_O2,
      bc_top = bc_top
    )


    # Recalculate redox state after oxygen diffusion
    anox <- anoxia_field(
      O2,
      o2_crit = o2_crit,
      sharpness = anoxia_sharpness
    )


    # --------------------------------------------------------
    # Biofilm development
    # --------------------------------------------------------

    B <- B +
      anox *
      k_bio *
      C *
      phi *
      (1 - B)

    B <- pmax(
      pmin(B, 1),
      0
    )


    # --------------------------------------------------------
    # Tissue degradation
    # --------------------------------------------------------
    #
    # Oxic contribution:
    #
    #   (1 - anox)
    #
    # Anaerobic contribution:
    #
    #   anaerobic_decay_ratio * anox
    #
    # This means complete anoxia no longer automatically
    # implies complete cessation of degradation.
    # --------------------------------------------------------

    redox_decay_multiplier <-
      (1 - anox) +
      anaerobic_decay_ratio * anox


    # Biofilm protection
    #
    # biofilm_protection = 0
    #   no protective effect
    #
    # biofilm_protection = 0.7
    #   maximum 70% reduction in degradation
    #
    biofilm_multiplier <-
      1 - biofilm_protection * B


    decay <-
      k_tissue *
      C *
      redox_decay_multiplier *
      biofilm_multiplier


    C <- C - decay

    C <- pmax(
      pmin(C, 1),
      0
    )


    # --------------------------------------------------------
    # Oxygen consumption
    # --------------------------------------------------------

    O2 <- O2 -
      k_resp_aero * C * (1 - anox) -
      k_resp_anaer * B * anox -
      k_sed


    O2 <- pmax(
      pmin(O2, 1),
      0
    )

    O2[1, ] <- bc_top


    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------

    if (record_history) {

      anox_diag <- anoxia_field(
        O2,
        o2_crit = o2_crit,
        sharpness = anoxia_sharpness
      )

      mean_anoxia[t] <-
        mean(anox_diag[body_mask])

      mean_o2[t] <-
        mean(O2[body_mask])

      mean_biofilm[t] <-
        mean(B[body_mask])

      total_survival[t] <-
        compute_tissue_survival(
          C,
          C_initial,
          body_mask
        )

      tissue_survival[t, ] <-
        compute_tissue_specific_survival(
          C,
          C_initial,
          tissue_ref
        )
    }
  }


  # ==========================================================
  # Final diagnostics
  # ==========================================================

  anox_final <- anoxia_field(
    O2,
    o2_crit = o2_crit,
    sharpness = anoxia_sharpness
  )


  final_survival <-
    compute_tissue_survival(
      C,
      C_initial,
      body_mask
    )


  final_tissue_survival <-
    compute_tissue_specific_survival(
      C,
      C_initial,
      tissue_ref
    )


  final_anoxia <-
    mean(anox_final[body_mask])


  final_o2 <-
    mean(O2[body_mask])


  final_biofilm <-
    mean(B[body_mask])


  # ----------------------------------------------------------
  # Cumulative anoxic exposure
  #
  # Mean anoxia across the body through simulation time.
  # Normalised to 0-1.
  # ----------------------------------------------------------

  cumulative_anoxia <-
    if (record_history) {
      mean(mean_anoxia)
    } else {
      NA_real_
    }


  # ----------------------------------------------------------
  # Time to predominantly anoxic conditions
  #
  # Defined provisionally as the first time mean body anoxia
  # exceeds 0.5.
  #
  # This is a diagnostic rather than a definition of
  # preservation and can therefore be changed during analysis.
  # ----------------------------------------------------------

  time_to_anoxia <-
    if (record_history &&
        any(mean_anoxia >= 0.5)) {

      which(mean_anoxia >= 0.5)[1]

    } else {

      NA_integer_
    }


  # ==========================================================
  # Return simulation
  # ==========================================================

  list(

    # Primary preservation outputs
    final_survival = final_survival,
    final_tissue_survival = final_tissue_survival,

    # Redox diagnostics
    final_anoxia = final_anoxia,
    final_o2 = final_o2,
    cumulative_anoxia = cumulative_anoxia,
    time_to_anoxia = time_to_anoxia,

    # Biofilm diagnostic
    final_biofilm = final_biofilm,

    # Final spatial fields
    C_final = C,
    O2_final = O2,
    anox_final = anox_final,
    B_final = B,

    # Initial state
    C_initial = C_initial,
    body_mask = body_mask,
    tissue_ref = tissue_ref,

    # Time series
    mean_anoxia = mean_anoxia,
    mean_o2 = mean_o2,
    mean_biofilm = mean_biofilm,
    total_survival = total_survival,
    tissue_survival = tissue_survival
  )
}
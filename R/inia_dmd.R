# ------------------------------------------------------------------------------
# Density management diagram (DMD) plotting
#
# Reference lines:
#   - Reineke's maximum (RD = 1.0), SDImax = 1250 for E. grandis subtropical SA
#     (Rachid-Casnati et al. 2024).
#   - Drew & Flewelling (1979) self-thinning lower limit (RD = 0.60).
#   - Drew & Flewelling crown-closure / competition-onset limit (RD = 0.35).
# Reference lines are literature-sourced, independent of `simulate_inia()`.
# ------------------------------------------------------------------------------


#' Density management diagram for a simulated stand
#'
#' Plots a simulated INIA trajectory on a Reineke-style stand density
#' diagram (log N versus log Dq) with literature-sourced reference lines
#' for maximum density (`RD = 1.0`), the self-thinning lower limit
#' (`RD = 0.60`), and the onset of competition (`RD = 0.35`).
#'
#' The reference lines are independent of `simulate_inia()`: they come
#' from `SDImax = 1250` (Rachid-Casnati et al. 2024) and Reineke's
#' `beta = 1.605`. The trajectory is overlaid with age labels at the start,
#' end, and any thinning events. Thinning events appear as downward
#' arrows from pre-thin to post-thin density.
#'
#' @param sim_result Output of [simulate_inia()].
#' @param label_ages Optional integer vector of ages to label on the
#'   trajectory. If `NULL`, labels start age, end age, and thinning ages.
#' @param dq_range Numeric length-2 vector, x-axis range in cm. If `NULL`,
#'   inferred from the trajectory with padding.
#' @param n_range Numeric length-2 vector, y-axis range in TPH. If `NULL`,
#'   inferred from the trajectory and the reference lines.
#' @return A ggplot object.
#' @examples
#' sim <- simulate_inia(
#'   SI = 28, N0 = 900, G0 = 7.0,
#'   t0 = 3, t_end = 18, zone = 7,
#'   thins = list(list(age = 7, N_after = 500))
#' )
#' inia_dmd_plot(sim)
#' @export
inia_dmd_plot <- function(sim_result, label_ages = NULL,
                          dq_range = NULL, n_range = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`ggplot2` is required for inia_dmd_plot(). ",
         "Install with install.packages(\"ggplot2\").")
  }

  traj <- sim_result$trajectory
  thins <- sim_result$thinnings

  p <- .inia_params$sdi
  SDImax <- p$SDImax; beta <- p$beta; Dq_ref <- p$Dq_ref

  if (is.null(dq_range)) {
    dq_range <- range(traj$DAP_medio) * c(0.7, 1.3)
    dq_range[1] <- max(dq_range[1], 5)
  }
  dq_grid <- exp(seq(log(dq_range[1]), log(dq_range[2]), length.out = 80))

  ref_levels <- c(
    "Reineke maximum (RD = 1.0)"          = 1.00,
    "Self-thinning lower limit (RD = 0.60)" = 0.60,
    "Competition onset (RD = 0.35)"       = 0.35
  )
  ref <- do.call(rbind, lapply(seq_along(ref_levels), function(i) {
    rd <- ref_levels[i]
    data.frame(
      Dq    = dq_grid,
      N     = SDImax * rd * (Dq_ref / dq_grid)^beta,
      label = factor(names(ref_levels)[i], levels = names(ref_levels))
    )
  }))

  if (is.null(n_range)) {
    # Bound the y-axis by the data and the Reineke-max line evaluated AT
    # the trajectory's Dq range. This keeps the headroom-to-max visible
    # without leaving a large empty region above small-Dq sections of
    # the reference lines.
    ref_max_at_traj <- SDImax * (Dq_ref / range(traj$DAP_medio))^beta
    n_range <- range(c(traj$N, ref_max_at_traj)) * c(0.85, 1.15)
  }

  if (is.null(label_ages)) {
    label_ages <- unique(c(min(traj$age), max(traj$age),
                           if (!is.null(thins)) thins$age))
  }
  labs <- traj[traj$age %in% label_ages, c("age", "DAP_medio", "N")]

  thin_segs <- if (!is.null(thins)) {
    data.frame(
      Dq_pre  = thins$Dq_pre,  N_pre  = thins$N_pre,
      Dq_post = thins$Dq_post, N_post = thins$N_post
    )
  } else NULL

  gg <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = ref,
      ggplot2::aes(x = .data$Dq, y = .data$N, linetype = .data$label),
      colour = "grey40"
    ) +
    ggplot2::geom_path(
      data = traj,
      ggplot2::aes(x = .data$DAP_medio, y = .data$N),
      colour = "steelblue", linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = traj,
      ggplot2::aes(x = .data$DAP_medio, y = .data$N),
      colour = "steelblue", size = 1.6
    )

  if (!is.null(thin_segs)) {
    gg <- gg + ggplot2::geom_segment(
      data = thin_segs,
      ggplot2::aes(x = .data$Dq_pre, xend = .data$Dq_post,
                   y = .data$N_pre, yend = .data$N_post),
      colour = "firebrick", linewidth = 0.7,
      arrow = ggplot2::arrow(length = ggplot2::unit(2.5, "mm"), type = "closed")
    )
  }

  gg +
    ggplot2::geom_text(
      data = labs,
      ggplot2::aes(x = .data$DAP_medio, y = .data$N,
                   label = paste0("age ", .data$age)),
      nudge_y = 0.05, vjust = -0.5, size = 3, colour = "black"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(xlim = dq_range, ylim = n_range) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "dashed", "dotted"), name = NULL
    ) +
    ggplot2::labs(
      x = "Quadratic mean diameter (cm)",
      y = "Trees per hectare",
      title = "Density management diagram"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ------------------------------------------------------------------------------
# Head-to-head comparison harness for simulate_inia() and simulate_grandis()
#
# Runs both simulators on a single site spec and returns a tidy data frame
# of the two trajectories side-by-side. A small plotting helper overlays
# common state variables for quick visual inspection. Designed for the
# "what does each model predict at this site?" workflow; does NOT require
# PSP calibration data.
# ------------------------------------------------------------------------------


#' Compare `simulate_inia()` and `simulate_grandis()` at a single site
#'
#' Runs both simulators on shared stand inputs (SI, N0, G0, ages,
#' thinning regime) and the two simulator-specific site arguments
#' (`zone` for INIA, `PASW` / `Elev` / `slope` / `aspect` for the
#' augmented simulator), and returns the two trajectories aligned in
#' long format for direct inspection or plotting.
#'
#' @param SI Site index (m).
#' @param N0,G0 Initial trees/ha and basal area (m^2/ha).
#' @param zone INIA zone (7, 8, 9). Default 7.
#' @param PASW,Elev,slope,aspect Augmented-simulator site variables.
#'   See [simulate_grandis()] for definitions.
#' @param Hd0,dmax0,SDd0,t0,t_end,thins,dt Shared simulator inputs;
#'   passed to both simulators verbatim.
#' @param label Optional character label attached to every row of the
#'   returned `comparison` frame -- useful when stacking multiple
#'   sites with `do.call(rbind, ...)`.
#' @return A list with:
#' \describe{
#'   \item{`comparison`}{Long-format data frame with one row per
#'     (model, age) pair. Columns: `label` (if supplied), `model`
#'     ("INIA" or "augmented"), plus every state column from the
#'     simulator trajectories.}
#'   \item{`inia`}{Full result of [simulate_inia()].}
#'   \item{`grandis`}{Full result of [simulate_grandis()].}
#' }
#' @examples
#' cmp <- compare_inia_grandis(
#'   SI = 28, N0 = 900, G0 = 7,
#'   zone = 7, PASW = 130, Elev = 130, slope = 0, aspect = 0,
#'   Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
#'   t0 = 2, t_end = 18,
#'   label = "Z7 baseline"
#' )
#' head(cmp$comparison)
#' @export
compare_inia_grandis <- function(SI, N0, G0,
                                 zone = 7,
                                 PASW = 130, Elev = 130,
                                 slope = 0, aspect = 0,
                                 Hd0 = NULL, dmax0 = NULL, SDd0 = NULL,
                                 t0 = 1, t_end = 16,
                                 thins = NULL,
                                 dt = 1,
                                 label = NULL) {

  inia <- simulate_inia(SI = SI, N0 = N0, G0 = G0,
                        Hd0 = Hd0, dmax0 = dmax0, SDd0 = SDd0,
                        t0 = t0, t_end = t_end, zone = zone,
                        thins = thins, dt = dt)

  grandis <- simulate_grandis(SI = SI, N0 = N0, G0 = G0,
                              PASW = PASW, Elev = Elev,
                              slope = slope, aspect = aspect,
                              Hd0 = Hd0, dmax0 = dmax0, SDd0 = SDd0,
                              t0 = t0, t_end = t_end,
                              thins = thins, dt = dt)

  cols <- c("age", "AMD", "N", "AB", "DAP_medio", "DAP_max",
            "Desvio_DAP", "Vol_Total", "MAI", "SDI", "RD", "ICA")
  inia_df    <- inia$trajectory[, cols, drop = FALSE]
  grandis_df <- grandis$trajectory[, cols, drop = FALSE]

  inia_df$model    <- "INIA"
  grandis_df$model <- "augmented"

  comparison <- rbind(inia_df, grandis_df)
  if (!is.null(label)) {
    comparison$label <- label
    comparison <- comparison[, c("label", "model", cols)]
  } else {
    comparison <- comparison[, c("model", cols)]
  }
  rownames(comparison) <- NULL

  list(
    comparison = comparison,
    inia       = inia,
    grandis    = grandis
  )
}


#' Quick overlay plot of INIA vs the augmented simulator
#'
#' Plots `simulate_inia()` and `simulate_grandis()` trajectories on the
#' same axes for the variables in `vars`, faceting if `label` is
#' present in the input frame. Returns a ggplot.
#'
#' @param comparison Long-format frame returned in the
#'   `compare_inia_grandis()` `$comparison` element, or several such
#'   frames row-bound together (in which case the `label` column
#'   drives the facets).
#' @param vars Character vector of state variables to plot. Default
#'   `c("AMD", "AB", "Vol_Total", "N")`.
#' @return A ggplot object.
#' @examples
#' cmp <- compare_inia_grandis(
#'   SI = 28, N0 = 900, G0 = 7,
#'   zone = 7, PASW = 130, Elev = 130,
#'   Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
#'   t0 = 2, t_end = 18, label = "site"
#' )
#' plot_inia_grandis_compare(cmp$comparison)
#' @export
plot_inia_grandis_compare <- function(comparison,
                                       vars = c("AMD", "AB",
                                                "Vol_Total", "N")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`ggplot2` is required for plot_inia_grandis_compare(). ",
         "Install with install.packages(\"ggplot2\").")
  }

  pretty_labels <- c(
    AMD = "Dominant height (m)",
    AB = "Basal area (m^2/ha)",
    Vol_Total = "Total volume (m^3/ha)",
    N = "Trees per hectare",
    DAP_medio = "Quadratic mean diameter (cm)",
    DAP_max = "Maximum diameter (cm)",
    Desvio_DAP = "Diameter SD (cm)",
    MAI = "MAI (m^3/ha/yr)",
    SDI = "Stand density index",
    RD = "Relative density",
    ICA = "Annual increment (m^3/ha/yr)"
  )

  keep <- c(intersect(c("label", "model", "age"), names(comparison)), vars)
  df <- comparison[, keep, drop = FALSE]

  # Long form for facetting by variable
  long <- do.call(rbind, lapply(vars, function(v) {
    out <- df[, c(intersect(c("label", "model", "age"), names(df)), v),
              drop = FALSE]
    out$variable <- pretty_labels[v]
    if (is.na(out$variable[1])) out$variable <- v
    names(out)[ncol(out) - 1] <- "value"
    out
  }))
  long$variable <- factor(long$variable,
                          levels = unname(pretty_labels[vars]))

  facet_formula <- if ("label" %in% names(long)) {
    ggplot2::vars(.data$label, .data$variable)
  } else {
    ggplot2::vars(.data$variable)
  }

  ggplot2::ggplot(long,
                  ggplot2::aes(x = .data$age, y = .data$value,
                               colour = .data$model,
                               linetype = .data$model)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(facet_formula, scales = "free_y") +
    ggplot2::scale_colour_manual(
      values = c(INIA = "firebrick", augmented = "steelblue"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c(INIA = "dashed", augmented = "solid"),
      name = NULL
    ) +
    ggplot2::labs(x = "Age (years)", y = NULL,
                  title = "INIA vs augmented simulator") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
}

# ------------------------------------------------------------------------------
# Aboveground biomass, carbon, and CO2-equivalent estimation.
#
# Tree-level allometry integrated across the recovered inverse Weibull
# diameter distribution. Calibrated against SAG grandis 2021 Biomasa
# output (3 reference scenarios, 15 data points; RRMSE 5.4%). See
# Hirigoyen et al. 2021 (BOSQUE 42(1): 53-66) for the additive biomass
# system that motivated the allometric form, and Resquin et al. 2019
# (For Ecol Manag 438: 63-74) for E. grandis density / wood-density
# context. Carbon and CO2 conversions follow IPCC 2006 defaults for
# subtropical hardwoods (carbon fraction 0.49).
# ------------------------------------------------------------------------------


#' Individual-tree total height from DBH and dominant height
#'
#' Exponential saturation curve anchored so `h -> Hd` as `d -> dmax` and
#' `h -> 1.3` as `d -> 0`:
#'
#' \deqn{h(d) = 1.3 + (Hd - 1.3) \cdot
#'   \frac{1 - e^{-k d}}{1 - e^{-k \, d_{max}}}}
#'
#' @param d DBH (cm), scalar or vector.
#' @param Hd Dominant height of the stand (m).
#' @param dmax Maximum diameter in the stand (cm).
#' @param k Curvature parameter. Default `0.1482` from joint fit to the
#'   SAG 2021 biomass scenarios.
#' @return Total height (m). Same length as `d`. Floored at 1.3 m.
#' @family biomass
#' @seealso [inia_tree_agb()], [inia_stand_agb()], [inia_diam_dist()]
#' @examples
#' inia_tree_height(d = c(10, 25, 40), Hd = 28, dmax = 49)
#' @export
inia_tree_height <- function(d, Hd, dmax,
                             k = .inia_params$biomass$hd_k) {
  h <- 1.3 + (Hd - 1.3) * (1 - exp(-k * d)) / (1 - exp(-k * dmax))
  pmax(h, 1.3)
}


#' Individual-tree aboveground biomass (kg dry matter)
#'
#' Allometric equation \eqn{AGB = a \cdot d^{b_1} \cdot h^{b_2}}
#' (kg per tree, total aboveground including stem, branches, and foliage).
#' Coefficients (a, b1, b2) jointly fitted with the height-diameter
#' curvature parameter against three SAG 2021 scenarios; see package
#' calibration notes.
#'
#' @param d DBH (cm), scalar or vector.
#' @param h Total height (m), same length as `d`.
#' @return AGB (kg per tree).
#' @family biomass
#' @seealso [inia_tree_height()], [inia_stand_agb()]
#' @examples
#' inia_tree_agb(d = 35, h = 28)
#' inia_tree_agb(d = c(20, 30, 40), h = c(22, 27, 31))
#' @export
inia_tree_agb <- function(d, h) {
  p <- .inia_params$biomass
  p$tree_a * d^p$tree_b1 * h^p$tree_b2
}


#' Stand-level aboveground biomass, carbon, and CO2-equivalent
#'
#' Recovers the inverse Weibull diameter distribution from the supplied
#' stand-level moments, estimates per-class height with
#' [inia_tree_height()] and per-class biomass with [inia_tree_agb()],
#' then sums to the stand level.
#'
#' @param N Trees per hectare.
#' @param Dq Quadratic mean diameter (cm).
#' @param dmax Maximum diameter (cm).
#' @param SDd Standard deviation of diameters (cm).
#' @param Hd Dominant height (m).
#' @return Named list rounded to 1 decimal place:
#' \describe{
#'   \item{`agb`}{Aboveground dry biomass (t/ha).}
#'   \item{`carbon`}{Carbon (t C/ha) = `agb` * 0.49.}
#'   \item{`co2eq`}{CO2-equivalent (t CO2/ha) = `agb` * 1.797.}
#' }
#' @family biomass
#' @seealso [inia_add_biomass()] for trajectory-wide application.
#' @examples
#' # Reference scenario at age 10 (Zone 7, SI=30, 550 TPH)
#' inia_stand_agb(N = 364, Dq = 38.3, dmax = 49.2, SDd = 6.1, Hd = 30.0)
#' @export
inia_stand_agb <- function(N, Dq, dmax, SDd, Hd) {
  dd <- inia_diam_dist(N = N, Dq = Dq, dmax = dmax, SDd = SDd)
  if (nrow(dd) == 0) {
    return(list(agb = 0, carbon = 0, co2eq = 0))
  }

  h <- inia_tree_height(dd$class_mid, Hd = Hd, dmax = dmax)
  agb_kg <- inia_tree_agb(dd$class_mid, h)
  agb_t_ha <- sum(dd$freq * agb_kg) / 1000

  p <- .inia_params$biomass
  list(
    agb    = round(agb_t_ha, 1),
    carbon = round(agb_t_ha * p$carbon_fraction, 1),
    co2eq  = round(agb_t_ha * p$co2_factor,      1)
  )
}


#' Append biomass / carbon / CO2 columns to an INIA simulation
#'
#' Convenience wrapper that walks the trajectory of a [simulate_inia()]
#' result and adds `Biomasa` (t/ha), `Carbon` (t C/ha), and `CO2eq`
#' (t CO2/ha) columns by calling [inia_stand_agb()] on each row.
#'
#' @param sim Result of [simulate_inia()].
#' @return The same `sim` object with three columns appended to
#'   `sim$trajectory`. All other elements are unchanged.
#' @family biomass
#' @seealso [inia_stand_agb()], [simulate_inia()]
#' @examples
#' sim <- simulate_inia(
#'   SI = 30, N0 = 550, G0 = 1.7,
#'   Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
#'   t0 = 1, t_end = 16, zone = 7
#' )
#' sim <- inia_add_biomass(sim)
#' head(sim$trajectory[, c("age", "Vol_Total", "Biomasa", "Carbon", "CO2eq")])
#' @export
inia_add_biomass <- function(sim) {
  traj <- sim$trajectory
  rows <- lapply(seq_len(nrow(traj)), function(i) {
    inia_stand_agb(
      N    = traj$N[i],
      Dq   = traj$DAP_medio[i],
      dmax = traj$DAP_max[i],
      SDd  = traj$Desvio_DAP[i],
      Hd   = traj$AMD[i]
    )
  })
  traj$Biomasa <- vapply(rows, `[[`, numeric(1), "agb")
  traj$Carbon  <- vapply(rows, `[[`, numeric(1), "carbon")
  traj$CO2eq   <- vapply(rows, `[[`, numeric(1), "co2eq")
  sim$trajectory <- traj
  sim
}

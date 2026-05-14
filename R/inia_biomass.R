# ------------------------------------------------------------------------------
# Aboveground biomass, carbon, and CO2-equivalent estimation.
#
# Tree-level allometric equations from Winck et al. 2015 (Ciencia
# Florestal 25(3): 595-606) for E. grandis in northeastern Argentina
# (Misiones + Corrientes; n=41 destructively sampled trees, ages 4-32 yr,
# DBH 16.5-38.1 cm, h 18.4-51.0 m). The default is the Zone 1 (Misiones
# Norte/Alta) fit, which has the broader age range and larger sample
# (n=23, ages 4-32). Stand-level biomass is obtained by integrating the
# per-tree equations across the recovered inverse Weibull diameter
# distribution. Carbon and CO2-equivalent conversions follow IPCC 2006
# defaults for subtropical hardwoods (carbon fraction 0.49).
#
# Audit history (see CLAUDE.md "Biomass" notes):
#   - v0.1.0  : single equation calibrated against SAG 2021 Biomasa
#               output. Realistic densities but not literature-traceable.
#   - v0.2.0  : Winck et al. 2015 NE Argentina E. grandis (current).
#               Within <5% of v0.1.0 predictions across the design
#               space, with full literature provenance and physical
#               wood-density consistency (~425-555 kg/m^3).
#               Hirigoyen et al. 2021 BOSQUE was evaluated but rejected
#               on the basis of an implausibly high implied stem
#               density (~1100 kg/m^3, roughly 2x physical wood density).
# ------------------------------------------------------------------------------


# Internal helper: evaluate B = fc * exp(b0 + b1*ln(d) + b2*ln(h)).
.winck_eqn <- function(d, h, p) {
  p$fc * exp(p$b0 + p$b1 * log(d) + p$b2 * log(h))
}


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
#' @param k Curvature parameter. Default `0.1482`. The h-d curve is an
#'   internal helper for estimating per-class h from a recovered diameter
#'   distribution; it is operationally independent of the Winck 2015
#'   biomass allometry, which takes h directly as an input.
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


#' Individual-tree stem biomass (kg dry matter)
#'
#' Log-linear allometric equation
#' \eqn{B_{stem} = f_c \cdot \exp(b_0 + b_1 \ln d + b_2 \ln h)}
#' for *E. grandis* stem (over-bark) dry biomass, kg per tree, from
#' Winck et al. (2015) Table 3 (Zone 1, Misiones Norte/Alta, Argentina;
#' n=23 destructively sampled trees, ages 4-32 yr, DBH 16.5-38.1 cm,
#' h 18.4-51.0 m). R^2 = 0.99, ECMP = 0.02; `fc` = 1.01 is the Meyer
#' (1941) correction for log-bias.
#'
#' @param d DBH (cm), scalar or vector.
#' @param h Total height (m), same length as `d`.
#' @return Stem dry biomass (kg per tree).
#' @references Winck RA, Fassola HE, Barth SR, Crechi EH, Keller AE,
#'   Videla D, Zaderenko C (2015). Modelos predictivos de biomasa aerea
#'   de *Eucalyptus grandis* para el noreste de Argentina. *Ciencia
#'   Florestal* 25(3): 595-606.
#' @family biomass
#' @examples
#' inia_tree_stem(d = 30, h = 28)
#' @export
inia_tree_stem <- function(d, h) {
  .winck_eqn(d, h, .inia_params$biomass$stem)
}


#' Individual-tree total-branches biomass (kg dry matter)
#'
#' Log-linear allometric equation
#' \eqn{B_{branches} = f_c \cdot \exp(b_0 + b_1 \ln d + b_2 \ln h)}
#' for *E. grandis* total branches (small + large) dry biomass, kg per
#' tree, from Winck et al. (2015) Table 3 (Zone 1, NE Argentina).
#' R^2 = 0.83, ECMP = 0.25; `fc` = 1.13 is the Meyer log-bias correction.
#'
#' @param d DBH (cm), scalar or vector.
#' @param h Total height (m), same length as `d`.
#' @return Branches dry biomass (kg per tree).
#' @references Winck RA et al. (2015). *Ciencia Florestal* 25(3): 595-606.
#' @family biomass
#' @examples
#' inia_tree_branches(d = 30, h = 28)
#' @export
inia_tree_branches <- function(d, h) {
  .winck_eqn(d, h, .inia_params$biomass$branches)
}


#' Individual-tree aboveground biomass (kg dry matter)
#'
#' Log-linear allometric equation
#' \eqn{AGB = f_c \cdot \exp(b_0 + b_1 \ln d + b_2 \ln h)}
#' for *E. grandis* total aboveground dry biomass (stem + branches +
#' foliage), kg per tree, from Winck et al. (2015) Table 3 (Zone 1, NE
#' Argentina). R^2 = 0.99, ECMP = 0.02; `fc` = 1.01.
#'
#' The total-AGB equation is fitted independently of the stem and
#' branches equations in the source paper; the implied component
#' breakdown (stem + branches predicted separately) differs from this
#' total by a small residual (~2-4%) that corresponds to leaves and
#' fit noise. Use [inia_tree_stem()] and [inia_tree_branches()] for
#' component-level estimates.
#'
#' @param d DBH (cm), scalar or vector.
#' @param h Total height (m), same length as `d`.
#' @return AGB (kg per tree).
#' @references Winck RA et al. (2015). *Ciencia Florestal* 25(3): 595-606.
#' @family biomass
#' @seealso [inia_tree_stem()], [inia_tree_branches()], [inia_stand_agb()]
#' @examples
#' inia_tree_agb(d = 35, h = 28)
#' inia_tree_agb(d = c(20, 30, 40), h = c(22, 27, 31))
#' @export
inia_tree_agb <- function(d, h) {
  .winck_eqn(d, h, .inia_params$biomass$agb)
}


#' Stand-level aboveground biomass, carbon, and CO2-equivalent
#'
#' Recovers the inverse Weibull diameter distribution from the supplied
#' stand-level moments, estimates per-class height with
#' [inia_tree_height()], computes per-class total / stem / branches
#' biomass via the Winck et al. 2015 equations, and sums to the stand
#' level.
#'
#' @param N Trees per hectare.
#' @param Dq Quadratic mean diameter (cm).
#' @param dmax Maximum diameter (cm).
#' @param SDd Standard deviation of diameters (cm).
#' @param Hd Dominant height (m).
#' @return Named list rounded to 1 decimal place:
#' \describe{
#'   \item{`agb`}{Total aboveground dry biomass (t/ha), from the
#'     independently-fitted total-AGB equation.}
#'   \item{`stem`}{Stem dry biomass (t/ha).}
#'   \item{`branches`}{Branches dry biomass (t/ha).}
#'   \item{`carbon`}{Carbon (t C/ha) = `agb` * 0.49.}
#'   \item{`co2eq`}{CO2-equivalent (t CO2/ha) = `agb` * 1.797.}
#' }
#' Note: `stem + branches` does not necessarily equal `agb` exactly
#' because the source-paper equations were fit independently rather
#' than under an additivity constraint. The residual (~2-4%) corresponds
#' to leaves and fit noise.
#' @family biomass
#' @seealso [inia_add_biomass()] for trajectory-wide application.
#' @examples
#' # Stand at mid-rotation
#' inia_stand_agb(N = 700, Dq = 28.0, dmax = 35.0, SDd = 4.5, Hd = 24.0)
#' @export
inia_stand_agb <- function(N, Dq, dmax, SDd, Hd) {
  dd <- inia_diam_dist(N = N, Dq = Dq, dmax = dmax, SDd = SDd)
  if (nrow(dd) == 0) {
    return(list(agb = 0, stem = 0, branches = 0, carbon = 0, co2eq = 0))
  }

  h        <- inia_tree_height(dd$class_mid, Hd = Hd, dmax = dmax)
  agb_kg   <- inia_tree_agb(dd$class_mid, h)
  stem_kg  <- inia_tree_stem(dd$class_mid, h)
  br_kg    <- inia_tree_branches(dd$class_mid, h)

  agb_t  <- sum(dd$freq * agb_kg)  / 1000
  stem_t <- sum(dd$freq * stem_kg) / 1000
  br_t   <- sum(dd$freq * br_kg)   / 1000

  p <- .inia_params$biomass
  list(
    agb      = round(agb_t,  1),
    stem     = round(stem_t, 1),
    branches = round(br_t,   1),
    carbon   = round(agb_t * p$carbon_fraction, 1),
    co2eq    = round(agb_t * p$co2_factor,      1)
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
#'   SI = 28, N0 = 900, G0 = 7.0,
#'   Hd0 = 7.0, dmax0 = 13.0, SDd0 = 1.8,
#'   t0 = 2, t_end = 16, zone = 7
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

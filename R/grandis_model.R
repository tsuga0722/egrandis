# ------------------------------------------------------------------------------
# Augmented simulator submodels for E. grandis (`simulate_grandis()`)
#
# Literature-faithful alternative to the SAG-2021-refitted forms in
# inia_model.R, intended for users who want explicit control over site
# variables (PASW, elevation, slope, aspect) and density-driven self-
# thinning. End-to-end provenance:
#
#   Hd      RC2019 Eqn 11 (Tab S4 row 11)
#   G       RC2019 Eqn 12 (Tab S5 row 12)
#   dmax    RC2019 Eqn 13 (Tab S6 row 13)
#   SDd     RC2019 Eqn 14 (Tab S7 row 14)
#   N (bg)  Methol 2003 Eqn 6 (Clutter-Jones, SI-dependent a)
#   N (st)  Reineke self-thinning ceiling (soft-logistic blend at RD = 0.60)
#
# Coefficients are taken verbatim from the RC2019 supplementary tables.
# The predictor mapping for each augmented equation was determined from
# the paper text (page 240: "PASW was significant for all state variables
# except for SDd, whereas aspect modified by slope was not significant
# for G"; page 242: "For dmax, site variables PASW and alpha_s included
# in the E. grandis' equation"), cross-checked against the Fig 5
# magnitudes.
#
# Convention for site variables:
#   PASW       potentially available soil water (mm)
#   Elev       elevation (m)
#   slope      slope (percent; flat ground = 0)
#   aspect     azimuth angle (radians from north; 0=N, pi/2=E, pi=S, 3pi/2=W)
#   alpha_s    slope * sin(aspect) -- East-West component, positive East
#   alpha_c    slope * cos(aspect) -- North-South component, positive North
#
# References:
#   Rachid-Casnati C, Mason E, Woollons R (2019). iForest 12: 237-245.
#   Methol R (2003). INIA Serie Tecnica 131, Eqn 6.
#   Reineke LH (1933). J Agric Res 46: 627-638.
#   Drew TJ, Flewelling JW (1979). For Sci 25: 518-532.
# ------------------------------------------------------------------------------


# Parameter table -------------------------------------------------------------
# Internal: not exported. Edit here to retune any submodel of the
# `simulate_grandis()` simulator. INIA parameters live in `.inia_params`
# (R/inia_model.R) and are not affected.
.grandis_params <- list(

  # Hd: RC2019 Eqn 11 (Tab S4 row 11). Augmented Johnson-Schumacher.
  # A = a0 + a1*PASW + a2*alpha_c + a3*alpha_s
  hd = list(a0 = 3.6014737,  a1 = 0.0029913,  a2 = 0.0087987,
            a3 = 0.0175036,  k  = 2.4973491),

  # G: RC2019 Eqn 12 (Tab S5 row 12). Augmented Schumacher(2) + thinning.
  # A = a0 + a1*PASW + a2*Elev
  # C = c0 + c1*(Na/Nb)/tt  (same thinning modifier as INIA base form)
  g  = list(a0 = 2.611299,   a1 = 0.0115832,  a2 = -0.0033998,
            c0 = 1.1164274,  c1 = -1.0353749),

  # dmax: RC2019 Eqn 13 (Tab S6 row 13). Augmented Schumacher.
  # A = a0 + a1*PASW + a2*alpha_s
  # C = c0 + c1*alpha_s
  dmax = list(a0 = 3.4624929,  a1 = 0.005749,   a2 = 0.0714808,
              c0 = 0.5225331,  c1 = -0.0101125),

  # SDd: RC2019 Eqn 14 (Tab S7 row 14). Augmented von Bertalanffy-Richards.
  # A = a0 + a1*alpha_s
  sdd = list(a0 = 8.863645,   a1 = 0.450453,   b = 0.069762),

  # Exogenous background mortality: a constant annual rate that captures
  # density-INDEPENDENT losses (wind, frost, drought, pests, lightning).
  # Competition mortality is added on top via the Reineke ceiling below.
  # Default 0.5%/year aligns with the residual mortality observed in
  # well-managed E. grandis plantations once the density-driven component
  # is accounted for (see e.g. Rachid-Casnati et al. 2024, which reports
  # 16% mortality over 18.6 yr at 810 TPH unthinned, of which roughly
  # half is density-driven).
  # Methol 2003 Eqn 6 (Clutter-Jones with SI-dependent a) was evaluated
  # as the background term but rejected: it double-counts competition
  # because its 2003 calibration data already contained density effects,
  # so layering Reineke on top produced implausible mortality (~78% in
  # 11 years at high RD).
  mortality = list(exo_rate = 0.005),

  # Reineke self-thinning ceiling. Drew & Flewelling (1979) lower limit
  # of self-thinning at RD = 0.60. SDImax = 1250 from Rachid-Casnati
  # et al. 2024 for E. grandis subtropical South America. The soft-
  # logistic blend uses `mort_k = 12` for a transition width of about
  # +/- 0.05 RD units around RD_50.
  reineke = list(SDImax = 1250, beta = 1.605, Dq_ref = 25,
                 RD_50 = 0.60, mort_k = 12),

  # Volume convention. Fang taper integration is over-bark by default;
  # apply `bark_factor` to convert to under-bark for direct comparison
  # with INIA Vol_Total. 0.82 is a published mean for E. grandis but
  # varies with tree size; exposed as a `simulate_grandis()` argument.
  volume = list(bark_factor = 0.82)
)


# Submodels -------------------------------------------------------------------
# All submodels are internal. Users access these through `simulate_grandis()`
# in R/grandis_simulate.R (Phase 2, forthcoming).

# Derive the East-West and North-South slope components from slope (in
# percent) and aspect (in radians from north). Returns a length-2 numeric
# vector c(alpha_s, alpha_c). Used by Hd, dmax, and SDd.
.grandis_slope_components <- function(slope = 0, aspect = 0) {
  c(alpha_s = slope * sin(aspect),
    alpha_c = slope * cos(aspect))
}


# Project dominant height from t1 to t2 under augmented Eqn 11.
# Hd1 in m, t1/t2 in years, PASW in mm, slope in %, aspect in radians.
grandis_hd <- function(Hd1, t1, t2, PASW, slope = 0, aspect = 0) {
  p <- .grandis_params$hd
  ac <- .grandis_slope_components(slope, aspect)
  A <- p$a0 + p$a1 * PASW + p$a2 * ac["alpha_c"] + p$a3 * ac["alpha_s"]
  R <- (t1 + p$k) / (t2 + p$k)
  unname(exp(log(Hd1) * R + A * (1 - R)))
}


# Convenience: dominant height at age t given site index SI at base_age.
grandis_hd_from_si <- function(t, SI, base_age = 10,
                               PASW, slope = 0, aspect = 0) {
  grandis_hd(SI, base_age, t, PASW = PASW, slope = slope, aspect = aspect)
}


#' Estimate site index from a dominant-height measurement (augmented form)
#'
#' Given a dominant-height observation `Hd` at age `age`, returns the
#' implied site index (dominant height at `base_age`, default 10 yr)
#' under the augmented RC2019 Eqn 11 height curve, parameterised by
#' the site variables `PASW`, `slope`, and `aspect`. Like
#' [inia_si_from_hd()], the calculation is exact because the
#' polymorphic Johnson-Schumacher form is path-independent.
#'
#' Vectorised over `Hd` and `age` (must be the same length); the site
#' variables are scalar (one site at a time). Loop over plots from
#' different sites by Map/lapply.
#'
#' Note on direction of the site-variable effect: the augmented Eqn 11
#' polymorphic curve is shifted upward at high PASW or favourable
#' aspect. For a measurement *before* `base_age` (e.g. age 8, base
#' age 10), recovering SI projects forward and the shifted curve
#' yields a *higher* SI estimate for the same (Hd, age) pair. For a
#' measurement *past* `base_age`, the projection runs backward and
#' the sign flips. Either way, feeding the recovered SI back into
#' [simulate_grandis()] with the same site variables reproduces the
#' observation.
#'
#' @param Hd Measured dominant height (m). Numeric, length >= 1.
#' @param age Age at measurement (years). Numeric, same length as `Hd`.
#' @param PASW Potentially available soil water (mm). Scalar.
#' @param slope Slope of the site (percent). Default 0 (flat).
#' @param aspect Aspect (azimuth) of the slope, radians from north
#'   (0 = N, pi/2 = E, pi = S, 3*pi/2 = W). Default 0.
#' @param base_age Reference age for site index (years). Default 10.
#' @return Estimated site index (m), same length as `Hd`.
#' @examples
#' # Single plot, with site variables
#' grandis_si_from_hd(Hd = 26, age = 8, PASW = 140, slope = 5, aspect = pi/4)
#'
#' # Multiple plots on the same site
#' grandis_si_from_hd(Hd  = c(22, 26, 30),
#'                    age = c( 7,  8, 10),
#'                    PASW = 140, slope = 5, aspect = pi/4)
#' @export
grandis_si_from_hd <- function(Hd, age,
                                PASW, slope = 0, aspect = 0,
                                base_age = 10) {
  if (length(Hd) != length(age)) {
    stop("`Hd` and `age` must have the same length")
  }
  if (any(Hd <= 0.1) || any(age <= 0)) {
    stop("`Hd` and `age` must be positive")
  }
  if (length(PASW) != 1 || length(slope) != 1 || length(aspect) != 1) {
    stop("`PASW`, `slope`, and `aspect` must be scalars (one site at a time)")
  }
  grandis_hd(Hd, t1 = age, t2 = base_age,
             PASW = PASW, slope = slope, aspect = aspect)
}


# Project basal area from t1 to t2 under augmented Eqn 12. The thinning
# modifier `c1 * (Na/Nb)/tt` matches the INIA base form: Na = trees/ha
# after thin, Nb = before, tt = age at thin.
grandis_g <- function(G1, t1, t2, PASW, Elev,
                      Na = NA, Nb = NA, tt = NA) {
  p <- .grandis_params$g
  A <- p$a0 + p$a1 * PASW + p$a2 * Elev
  C <- if (!is.na(Na) && !is.na(Nb) && !is.na(tt) && tt > 0) {
    p$c0 + p$c1 * (Na / Nb) / tt
  } else {
    p$c0
  }
  R <- (t1 / t2)^C
  exp(log(max(G1, 0.001)) * R + A * (1 - R))
}


# Project maximum diameter from t1 to t2 under augmented Eqn 13.
# Both A and the curvature exponent C are augmented with alpha_s.
grandis_dmax <- function(dmax1, t1, t2, PASW, slope = 0, aspect = 0) {
  p <- .grandis_params$dmax
  ac <- .grandis_slope_components(slope, aspect)
  A <- p$a0 + p$a1 * PASW + p$a2 * ac["alpha_s"]
  C <- p$c0 + p$c1 * ac["alpha_s"]
  R <- (t1 / t2)^C
  unname(exp(log(max(dmax1, 0.1)) * R + A * (1 - R)))
}


# Project standard deviation of diameters from t1 to t2 under
# augmented Eqn 14.
grandis_sdd <- function(SDd1, t1, t2, slope = 0, aspect = 0) {
  p <- .grandis_params$sdd
  ac <- .grandis_slope_components(slope, aspect)
  A <- p$a0 + p$a1 * ac["alpha_s"]
  ratio <- log(1 - exp(-p$b * t2)) / log(1 - exp(-p$b * t1))
  unname(A * (SDd1 / A)^ratio)
}


# Background mortality projection: constant annual exogenous-mortality
# rate, applied as continuous exponential decay. Returns the projected
# N from t1 to t2 BEFORE any Reineke self-thinning adjustment.
#
# `rate` is the fraction of trees lost to exogenous causes per year
# (default 0.005 = 0.5%/year, taken from .grandis_params$mortality).
# The continuous formulation makes the function time-step independent.
grandis_n_background <- function(N1, t1, t2, rate = NULL) {
  if (is.null(rate)) rate <- .grandis_params$mortality$exo_rate
  if (N1 <= 0) return(0)
  N1 * exp(-rate * (t2 - t1))
}


# Combined mortality: constant-rate exogenous background + Reineke
# self-thinning ceiling via soft logistic. The function applies the
# exogenous decay first, then derives the implied RD at the new Dq
# (computed from the projected G and the post-exogenous N), and
# blends N_bg toward the Reineke ceiling N_max(Dq) via a logistic in
# (RD_bg - RD_50). At low RD the function returns ~N_bg (pure
# background); at high RD it returns ~min(N_bg, N_ceiling).
grandis_n <- function(N1, t1, t2, G,
                      exo_rate = NULL,
                      RD_50    = NULL,
                      mort_k   = NULL,
                      SDImax   = NULL,
                      beta     = NULL,
                      Dq_ref   = NULL) {
  rp <- .grandis_params$reineke
  if (is.null(RD_50))   RD_50  <- rp$RD_50
  if (is.null(mort_k))  mort_k <- rp$mort_k
  if (is.null(SDImax))  SDImax <- rp$SDImax
  if (is.null(beta))    beta   <- rp$beta
  if (is.null(Dq_ref))  Dq_ref <- rp$Dq_ref

  N_bg <- grandis_n_background(N1, t1, t2, rate = exo_rate)
  if (N_bg <= 0 || G <= 0) return(N_bg)

  Dq_bg <- sqrt(G / N_bg * 40000 / pi)
  if (Dq_bg <= 0) return(N_bg)

  RD_bg     <- (N_bg * (Dq_bg / Dq_ref)^beta) / SDImax
  N_ceiling <- SDImax * RD_50 * (Dq_ref / Dq_bg)^beta

  weight <- 1 / (1 + exp(-mort_k * (RD_bg - RD_50)))
  N_at_ceil <- min(N_bg, N_ceiling)

  (1 - weight) * N_bg + weight * N_at_ceil
}


# Stand-level under-bark volume (m^3/ha) via Fang taper x recovered
# inverse-Weibull diameter distribution. Reuses INIA's diameter-recovery,
# h-d, and Fang total-volume helpers. Multiplies the over-bark integration
# by `bark_factor` (default 0.82 for E. grandis) to give an under-bark
# volume comparable to INIA's stand-level Vol_Total.
grandis_vol <- function(G, N, Hd, dmax, SDd, bark_factor = NULL) {
  if (is.null(bark_factor)) bark_factor <- .grandis_params$volume$bark_factor
  if (G <= 0.01 || N <= 0) return(0)

  Dq <- sqrt(G / N * 40000 / pi)
  if (Dq <= 0.1) return(0)

  dd <- inia_diam_dist(N = N, Dq = Dq, dmax = dmax, SDd = SDd)
  if (nrow(dd) == 0) return(0)

  h_class <- inia_height_class(dd$class_mid, Dq = Dq, Hd = Hd)
  v_per_tree <- vapply(seq_len(nrow(dd)), function(i) {
    if (h_class[i] <= 1.3 || dd$class_mid[i] <= 0) return(0)
    inia_tree_total_vol(D = dd$class_mid[i], H = h_class[i])
  }, numeric(1))

  sum(dd$freq * v_per_tree) * bark_factor
}

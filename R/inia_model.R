# ------------------------------------------------------------------------------
# INIA E. grandis 2021 - submodel functions and parameters
#
# Faithful reproduction of the SAG grandis 2021 online simulator
# (https://sag.inia.uy). Per-submodel validation status is summarised in
# the package README.
#
# References:
#   Rachid-Casnati C, Mason E, Woollons R (2019). iForest 12: 237-245.
#   Methol R (2003). INIA Serie Tecnica 131.
#   Methol R (2001). PhD thesis, University of Canterbury.
# ------------------------------------------------------------------------------


# Parameter table -------------------------------------------------------------
# Internal: not exported. Edit here to retune any submodel.
.inia_params <- list(

  # Hd: RC2019 eqn 5 (Johnson-Schumacher, polymorphic)
  hd = list(a0 = 4.00389, a1 = 0.19295, k = 3.06896),

  # G: RC2019 eqn 6 (Schumacher(2) + thinning modifier)
  g = list(a0 = 3.7534, a1 = 0.27345, c0 = 1.07956, c1 = -0.93323),

  # N: Clutter-Jones, refit to SAG 2021 output (zone-specific)
  n_z7  = list(a = -0.006319,  b = 0.003436,  c = 0.5669),
  n_z89 = list(a = -20.951726, b = 0.000000,  c = 5.4819),

  # V: stand-level exponential, refit to SAG 2021 output
  v = list(b0 = -1.1317, b1 = 1.0045, b2 = 0.1040, b3 = -0.0167),

  # Gpost: SAG 2003 eqn 9 (empirical thinning BA reduction)
  gpost = list(a = 1.1499, b = 0.9356, c = 1.5167, d = 0.9887),

  # dmax: Schumacher, refit to SAG 2021 output (only Zone 7 fitted)
  dmax_z7 = list(A = 4.760478, C = 0.495837),

  # SDd: RC2019 eqn 8 (von Bertalanffy-Richards)
  sdd = list(a0 = 9.230291, a1 = 2.424718, b = 0.054145),

  # Inverse Weibull polynomial coefficients (SAG 2003 eqn 18)
  weibull_kz = c(-0.22004032, -0.001433169, 0.150611381,
                 -0.078575996,  0.004305716, 0.008804944),

  # Aboveground biomass module. Calibrated against SAG grandis 2021
  # output (3 scenarios, 15 data points, RRMSE 5.4%). See R/inia_biomass.R.
  biomass = list(
    # Tree-level AGB allometry: AGB(kg) = a * d^b1 * h^b2
    tree_a  = 0.040403,
    tree_b1 = 2.5094,
    tree_b2 = 0.2362,
    # Height-diameter saturation curve parameter
    hd_k = 0.1482,
    # IPCC subtropical hardwood carbon fraction; CO2eq factor = 0.49 * 44/12.
    carbon_fraction = 0.49,
    co2_factor      = 0.49 * 44 / 12
  )
)


# Submodels -------------------------------------------------------------------
# All submodels are internal. Users access the model through simulate_inia().

# Project dominant height from t1 to t2.
# Hd1 in m, t1/t2 in years, Z7 = 1 for Zone 7 else 0.
inia_hd <- function(Hd1, t1, t2, Z7 = 1) {
  p <- .inia_params$hd
  A <- p$a0 + p$a1 * Z7
  R <- (t1 + p$k) / (t2 + p$k)
  exp(log(Hd1) * R + A * (1 - R))
}

# Convenience: dominant height at age t given site index SI at base_age.
inia_hd_from_si <- function(t, SI, base_age = 10, Z7 = 1) {
  inia_hd(SI, base_age, t, Z7)
}

# Project basal area from t1 to t2. Thinning history modifies C.
# Na = trees/ha after thin, Nb = before, tt = age at thin.
inia_g <- function(G1, t1, t2, Z7 = 1, Na = NA, Nb = NA, tt = NA) {
  p <- .inia_params$g
  A <- p$a0 + p$a1 * Z7
  C <- if (!is.na(Na) && !is.na(Nb) && !is.na(tt) && tt > 0) {
    p$c0 + p$c1 * (Na / Nb) / tt
  } else {
    p$c0
  }
  R <- (t1 / t2)^C
  exp(log(max(G1, 0.001)) * R + A * (1 - R))
}

# Project trees/ha (background mortality).
inia_n <- function(N1, t1, t2, Z7 = 1) {
  p <- if (Z7 == 1) .inia_params$n_z7 else .inia_params$n_z89
  val <- N1^p$a + p$b * ((t2 / 10)^p$c - (t1 / 10)^p$c)
  if (val <= 0) return(0)
  val^(1 / p$a)
}

# Stand-level under-bark volume (m3/ha).
inia_vol <- function(G, N, Hd, Z7 = 1) {
  p <- .inia_params$v
  Dg <- sqrt(G / max(N, 1) * 40000 / pi)
  if (Dg < 0.1 || G < 0.01) return(0)
  exp(p$b0 + p$b1 * log(G * Hd) + p$b2 * (Hd / Dg) + p$b3 * Z7)
}

# Post-thinning basal area (SAG 2003 eqn 9). Empirical.
inia_gpost <- function(Gant, Nant, Npost) {
  p <- .inia_params$gpost
  frac_removed <- 1 - Npost / Nant
  p$a * Gant^p$b * (1 - frac_removed^p$c)^p$d
}

# Project max diameter. Only Zone 7 calibrated.
inia_dmax <- function(dmax1, t1, t2, Z7 = 1) {
  p <- .inia_params$dmax_z7
  R <- (t1 / t2)^p$C
  exp(log(max(dmax1, 0.1)) * R + p$A * (1 - R))
}

# Project SD of diameters.
inia_sdd <- function(SDd1, t1, t2, Z7 = 1) {
  p <- .inia_params$sdd
  A <- p$a0 + p$a1 * Z7
  ratio <- log(1 - exp(-p$b * t2)) / log(1 - exp(-p$b * t1))
  A * (SDd1 / A)^ratio
}

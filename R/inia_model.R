# ------------------------------------------------------------------------------
# INIA E. grandis 2021 - submodel functions and parameters
#
# Faithful reproduction of the SAG grandis 2021 online simulator
# (https://sag.inia.uy). Per-submodel literature provenance audit (May 2026)
# is documented in CLAUDE.md "May 2026 deep-audit summary".
#
# Summary of submodel provenance after audit:
#   EXACT match to published source: Hd, G, SDd, Gpost, Weibull recovery
#   SAG-2021-refit (form from literature, coefficients re-derived): N, V, dmax
#
# References:
#   Rachid-Casnati C, Mason E, Woollons R (2019). iForest 12: 237-245.
#     Hd (Eqn 5 / Tab S4), G (Eqn 6 / Tab S5), dmax form (Eqn 7),
#     SDd (Eqn 8 / Tab S7).
#   Methol R (2003). INIA Serie Tecnica 131.
#     Gpost (Eqn 9), Weibull pdf and recovery (Eqns 10-11, 16-18).
#   Methol R (2001). PhD thesis, University of Canterbury.
#     CITASS grandis VBA implementation (different mortality and volume forms
#     than Methol 2003 documents -- the docs and code diverged even at that
#     vintage).
# ------------------------------------------------------------------------------


# Parameter table -------------------------------------------------------------
# Internal: not exported. Edit here to retune any submodel.
.inia_params <- list(

  # Hd: RC2019 eqn 5 (Johnson-Schumacher, polymorphic)
  hd = list(a0 = 4.00389, a1 = 0.19295, k = 3.06896),

  # G: RC2019 eqn 6 (Schumacher(2) + thinning modifier)
  g = list(a0 = 3.7534, a1 = 0.27345, c0 = 1.07956, c1 = -0.93323),

  # N: Clutter-Jones FORM (Methol 2003 Eqn 6), but coefficients REFIT to
  # SAG 2021 output -- the SAG 2021 system removed the site-dependent `a`
  # that Methol 2003 documented and also diverges from Methol's own 2001
  # CITASS VBA implementation. See CLAUDE.md "deep-audit summary".
  # Z7 max abs err ~12 trees vs SAG; Z8/9 ~20 trees (humped per-year-loss
  # shape that the monotonic Clutter-Jones form cannot fully fit).
  n_z7  = list(a = -0.006319,  b = 0.003436,  c = 0.5669),
  n_z89 = list(a = -0.199951,  b = 0.017897,  c = 0.955840),

  # V: stand-level Schumacher-Hall FORM (Methol 2003 Eqn 7), but
  # coefficients REFIT to SAG 2021 output. Note the b2 sign-flip between
  # Methol 2003 (b2 = -0.0761) and SAG 2021 refit (b2 = +0.104) -- this
  # is a real recalibration in SAG 2021, not an implementation bug.
  v = list(b0 = -1.1317, b1 = 1.0045, b2 = 0.1040, b3 = -0.0167),

  # Gpost: Methol 2003 Eqn 9 (empirical thinning BA reduction). EXACT match.
  gpost = list(a = 1.1499, b = 0.9356, c = 1.5167, d = 0.9887),

  # dmax: RC2019 Eqn 7 FORM (Schumacher), coefficients REFIT to SAG 2021.
  # RC2019 published values for Z7 (a0=4.269, a1=1.147, c0=0.454, c1=-0.159)
  # underestimate SAG 2021 dmax by 4-8 cm. Methol 2003 Eqn 14 used a
  # different (more complex) von Bertalanffy-Richards form with SI and N0
  # dependence; SAG 2021 simplified to the RC2019 form and re-fit.
  dmax_z7  = list(A = 4.760478, C = 0.495837),
  dmax_z89 = list(A = 4.055992, C = 0.650911),

  # SDd: RC2019 eqn 8 (von Bertalanffy-Richards)
  sdd = list(a0 = 9.230291, a1 = 2.424718, b = 0.054145),

  # Inverse Weibull polynomial coefficients (Methol 2003 Eqn 18). EXACT match.
  weibull_kz = c(-0.22004032, -0.001433169, 0.150611381,
                 -0.078575996,  0.004305716, 0.008804944),

  # Stand density index. Reineke's beta = 1.605 (1933 classical value);
  # SDImax = 1250 for E. grandis subtropical South America
  # (Rachid-Casnati et al. 2024).
  sdi = list(SDImax = 1250, beta = 1.605, Dq_ref = 25),

  # Aboveground biomass module. Tree-level prediction models from Winck
  # et al. 2015 (Ciencia Florestal 25(3): 595-606) for E. grandis in NE
  # Argentina (Misiones + Corrientes; n=41 destructively sampled trees,
  # ages 4-32 yr, DBH 16.5-38.1 cm, h 18.4-51.0 m). Zone 1 (Misiones N/A,
  # n=23, ages 4-32) is the default for its broader age range. Each
  # equation is of the form B = fc * exp(b0 + b1*ln(d) + b2*ln(h)),
  # where fc is the Meyer (1941) correction factor for the systematic
  # bias of logarithmic regressions. Components were fit independently
  # in the source paper (not via NSUR), so stem + branches != AGB
  # exactly -- the residual (~2-4%) corresponds to leaves and fit noise.
  # See R/inia_biomass.R.
  biomass = list(
    agb      = list(b0 = -3.36, b1 = 2.12, b2 =  0.65, fc = 1.01),
    stem     = list(b0 = -4.51, b1 = 1.83, b2 =  1.22, fc = 1.01),
    branches = list(b0 = -2.68, b1 = 3.73, b2 = -1.77, fc = 1.13),
    # Height-diameter saturation curve parameter (operationally separate
    # from the allometry above; used by inia_tree_height() to estimate
    # per-class h from a recovered diameter distribution).
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

# Project max diameter. Zone-specific.
inia_dmax <- function(dmax1, t1, t2, Z7 = 1) {
  p <- if (Z7 == 1) .inia_params$dmax_z7 else .inia_params$dmax_z89
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

# Reineke's stand density index. Independent of site/zone.
inia_sdi <- function(N, Dq) {
  p <- .inia_params$sdi
  N * (Dq / p$Dq_ref)^p$beta
}

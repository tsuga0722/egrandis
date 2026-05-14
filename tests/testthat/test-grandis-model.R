# Phase 1 tests for R/grandis_model.R
# Sensitivity / integration / Reineke-kick-in tests live in
# test-grandis-simulate.R (Phase 3, alongside the simulator).

# ---- helpers ----------------------------------------------------------------

test_that(".grandis_slope_components decomposes slope x aspect correctly", {
  # Flat ground -> both components zero regardless of aspect
  for (aspect in c(0, pi/2, pi, 3*pi/2)) {
    out <- egrandis:::.grandis_slope_components(slope = 0, aspect = aspect)
    expect_equal(unname(out["alpha_s"]), 0)
    expect_equal(unname(out["alpha_c"]), 0)
  }
  # 10% slope facing north -> alpha_c = +10, alpha_s = 0
  out <- egrandis:::.grandis_slope_components(slope = 10, aspect = 0)
  expect_equal(unname(out["alpha_s"]), 0,  tolerance = 1e-10)
  expect_equal(unname(out["alpha_c"]), 10, tolerance = 1e-10)
  # 10% slope facing east -> alpha_s = +10, alpha_c = 0
  out <- egrandis:::.grandis_slope_components(slope = 10, aspect = pi/2)
  expect_equal(unname(out["alpha_s"]), 10, tolerance = 1e-10)
  expect_equal(unname(out["alpha_c"]),  0, tolerance = 1e-10)
  # 10% slope facing west -> alpha_s = -10
  out <- egrandis:::.grandis_slope_components(slope = 10, aspect = 3*pi/2)
  expect_equal(unname(out["alpha_s"]), -10, tolerance = 1e-10)
})

# ---- Hd: RC2019 Eqn 11 ------------------------------------------------------

test_that("grandis_hd is an identity projection when t1 == t2", {
  expect_equal(grandis_hd(20, 10, 10, PASW = 120), 20)
})

test_that("grandis_hd returns plausible values across the data envelope", {
  # RC2019 E. grandis sites: hdom 4.4-46.7 m (Tab 2)
  for (PASW in c(80, 120, 160)) {
    Hd <- grandis_hd(7, 2, 16, PASW = PASW, slope = 3, aspect = 0)
    expect_gt(Hd, 15)   # any projection from age 2 to 16 produces >15m
    expect_lt(Hd, 60)   # ... and stays below E. grandis tall outlier
  }
})

test_that("grandis_hd responds correctly to PASW and aspect", {
  # Higher PASW -> taller Hd
  Hd_lo  <- grandis_hd(7, 2, 16, PASW = 80)
  Hd_hi  <- grandis_hd(7, 2, 16, PASW = 160)
  expect_gt(Hd_hi, Hd_lo)
  # NE-facing -> taller than SW-facing
  Hd_NE <- grandis_hd(7, 2, 16, PASW = 120, slope = 9, aspect = pi/4)
  Hd_SW <- grandis_hd(7, 2, 16, PASW = 120, slope = 9, aspect = 5*pi/4)
  expect_gt(Hd_NE, Hd_SW)
})

# ---- G: RC2019 Eqn 12 -------------------------------------------------------

test_that("grandis_g is an identity projection when t1 == t2", {
  expect_equal(grandis_g(20, 10, 10, PASW = 120, Elev = 150), 20,
               tolerance = 1e-10)
})

test_that("grandis_g responds correctly to PASW and Elev", {
  G_pasw_hi <- grandis_g(7, 2, 16, PASW = 160, Elev = 130)
  G_pasw_lo <- grandis_g(7, 2, 16, PASW =  80, Elev = 130)
  expect_gt(G_pasw_hi, G_pasw_lo)         # higher PASW -> higher G
  G_elev_lo <- grandis_g(7, 2, 16, PASW = 120, Elev =  50)
  G_elev_hi <- grandis_g(7, 2, 16, PASW = 120, Elev = 250)
  expect_gt(G_elev_lo, G_elev_hi)         # higher Elev -> lower G (a2 < 0)
})

test_that("grandis_g thinning modifier reduces post-thin curvature", {
  # With a thinning history, C = c0 + c1*(Na/Nb)/tt; c1 < 0 should slow
  # the approach toward the asymptote post-thin.
  G_unthinned <- grandis_g(20, 6, 12, PASW = 130, Elev = 100)
  G_thinned   <- grandis_g(20, 6, 12, PASW = 130, Elev = 100,
                           Na = 500, Nb = 900, tt = 5)
  expect_lt(G_thinned, G_unthinned)
})

# ---- dmax: RC2019 Eqn 13 ----------------------------------------------------

test_that("grandis_dmax responds to PASW and East-West aspect", {
  # Higher PASW -> larger dmax
  d_lo <- grandis_dmax(8, 2, 16, PASW = 116, slope = 5, aspect = pi/2)
  d_hi <- grandis_dmax(8, 2, 16, PASW = 169, slope = 5, aspect = pi/2)
  expect_gt(d_hi, d_lo)
  # East-facing -> larger dmax than West-facing
  d_E <- grandis_dmax(8, 2, 16, PASW = 140, slope = 9, aspect = pi/2)
  d_W <- grandis_dmax(8, 2, 16, PASW = 140, slope = 9, aspect = 3*pi/2)
  expect_gt(d_E, d_W)
})

# ---- SDd: RC2019 Eqn 14 -----------------------------------------------------

test_that("grandis_sdd responds to East-West aspect", {
  s_E <- grandis_sdd(1.3, 2, 16, slope = 5, aspect = pi/2)
  s_W <- grandis_sdd(1.3, 2, 16, slope = 5, aspect = 3*pi/2)
  expect_gt(s_E, s_W)
})

# ---- Mortality: Methol 2003 Eqn 6 + Reineke ceiling -------------------------

test_that("grandis_n_background gives sane mortality across the SI envelope", {
  # Calibration range: SI 22-35 (Methol 2003's a = 0.4577 - 0.0218*SI is
  # negative). Outside that range we clamp `a` to -0.005 to prevent the
  # documented anti-mortality artefact at low SI.
  for (SI in c(22, 25, 28, 30, 35)) {
    N16 <- grandis_n_background(N1 = 900, t1 = 2, t2 = 16, SI = SI)
    expect_lt(N16, 900, label = paste("mortality at SI =", SI))
    expect_gt(N16,   0, label = paste("non-negative N at SI =", SI))
  }
})

test_that("grandis_n_background clamps a at SI < ~22 to prevent anti-mortality", {
  # Below the clamp threshold, the formula must not produce N2 > N1.
  for (SI in c(15, 18, 20, 21)) {
    N16 <- grandis_n_background(N1 = 900, t1 = 2, t2 = 16, SI = SI)
    expect_lte(N16, 900, label = paste("no anti-mortality at SI =", SI))
  }
})

test_that("grandis_n collapses to background when RD is very low", {
  # Low N + small G -> Dq small -> RD tiny -> Reineke weight ~ 0 ->
  # function returns ~ N_background.
  N_bg <- grandis_n_background(N1 = 300, t1 = 2, t2 = 16, SI = 28)
  N_w  <- grandis_n(N1 = 300, t1 = 2, t2 = 16, SI = 28, G = 5)
  expect_equal(N_w, N_bg, tolerance = 0.5)
})

test_that("grandis_n pulls N down toward Reineke ceiling at high RD", {
  # High N + large G -> RD well above 0.6 -> heavy Reineke blend.
  # Note: because Dq depends on N and we use Dq_bg (pre-mortality) to
  # compute N_ceiling, the soft-logistic blend provides resistance
  # toward RD=0.60 but does not lock the post-mortality RD at that
  # value -- as N falls, Dq rises, partially restoring RD. The cap
  # below quantifies the resistance, not the equilibrium.
  N_bg <- grandis_n_background(N1 = 1500, t1 = 2, t2 = 16, SI = 28)
  N_w  <- grandis_n(N1 = 1500, t1 = 2, t2 = 16, SI = 28, G = 60)
  expect_lt(N_w, N_bg)
  # Substantial extra mortality is applied: > 25% of N_bg is removed.
  expect_lt(N_w, N_bg * 0.75)
  # The post-mortality RD sits between 0.60 and the pre-mortality value;
  # it does not exactly reach 0.60 because of the Dq-N feedback.
  Dq_w  <- sqrt(60 / N_w * 40000 / pi)
  Dq_bg <- sqrt(60 / N_bg * 40000 / pi)
  RD_w  <- (N_w  * (Dq_w  / 25)^1.605) / 1250
  RD_bg <- (N_bg * (Dq_bg / 25)^1.605) / 1250
  expect_gt(RD_w, 0.60)
  expect_lt(RD_w, RD_bg)
})

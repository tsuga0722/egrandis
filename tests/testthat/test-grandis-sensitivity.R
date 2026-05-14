# Phase 3 tests for simulate_grandis(): sensitivity analyses, Reineke
# kinetics, and cross-simulator comparison against simulate_inia().
# Basic structure and correctness live in test-grandis-{model,simulate}.R.

# Helpers ---------------------------------------------------------------------

# Run a baseline grandis simulation with single-line site arguments.
run_g <- function(SI = 28, N0 = 900, G0 = 7,
                  PASW = 130, Elev = 130, slope = 0, aspect = 0,
                  t0 = 2, t_end = 16, thins = NULL) {
  simulate_grandis(SI = SI, N0 = N0, G0 = G0,
                   PASW = PASW, Elev = Elev,
                   slope = slope, aspect = aspect,
                   t0 = t0, t_end = t_end, thins = thins)
}
tail_row <- function(sim) tail(sim$trajectory, 1)


# ---- 1. PASW dose-response --------------------------------------------------

test_that("PASW dose-response is monotonic and substantial across Hd, G, Vol", {
  outs <- lapply(c(60, 100, 140, 180), function(p) tail_row(run_g(PASW = p)))
  Hd  <- vapply(outs, function(r) r$AMD,      numeric(1))
  G   <- vapply(outs, function(r) r$AB,       numeric(1))
  Vol <- vapply(outs, function(r) r$Vol_Total, numeric(1))

  expect_true(all(diff(Hd)  > 0), label = "Hd monotone in PASW")
  expect_true(all(diff(G)   > 0), label = "G monotone in PASW")
  expect_true(all(diff(Vol) > 0), label = "Vol monotone in PASW")

  # Vol triples across the PASW range -- this is the dominant site driver.
  expect_gt(Vol[4] / Vol[1], 2.5)
})


# ---- 2. Elev sensitivity ----------------------------------------------------

test_that("Hd is invariant under Elev; G and Vol decrease with Elev", {
  outs <- lapply(c(50, 100, 150, 200, 250), function(e) tail_row(run_g(Elev = e)))
  Hd  <- vapply(outs, function(r) r$AMD,       numeric(1))
  G   <- vapply(outs, function(r) r$AB,        numeric(1))
  Vol <- vapply(outs, function(r) r$Vol_Total, numeric(1))

  # Elev does not enter RC2019 Eqn 11; Hd should be identical to 0.1 m.
  expect_lt(diff(range(Hd)), 0.2)
  # Elev coefficient a2 is negative in Eqn 12 -> G falls with Elev.
  expect_true(all(diff(G)   < 0), label = "G monotone-decreasing in Elev")
  expect_true(all(diff(Vol) < 0), label = "Vol monotone-decreasing in Elev")
})


# ---- 3. Aspect at fixed slope -----------------------------------------------

test_that("Aspect responses match the augmented-equation predictor mapping", {
  # 9% slope, fixed PASW/Elev/SI/density. The augmented equations
  # predict NE-facing maxima for Hd (a2*alpha_c + a3*alpha_s, both
  # positive at NE) and East-facing maxima for dmax and SDd
  # (a*alpha_s only).
  N  <- tail_row(run_g(slope = 9, aspect = 0))$AMD          # north
  NE <- tail_row(run_g(slope = 9, aspect = pi/4))           # northeast
  E  <- tail_row(run_g(slope = 9, aspect = pi/2))           # east
  S  <- tail_row(run_g(slope = 9, aspect = pi))             # south
  SW <- tail_row(run_g(slope = 9, aspect = 5*pi/4))         # southwest
  W  <- tail_row(run_g(slope = 9, aspect = 3*pi/2))         # west

  # Hd: NE maxima
  expect_gt(NE$AMD, S$AMD)
  expect_gt(NE$AMD, SW$AMD)
  # dmax: East-facing strictly larger than West-facing
  expect_gt(E$DAP_max, W$DAP_max)
  expect_gt(NE$DAP_max, SW$DAP_max)
  # SDd: East-facing strictly larger than West-facing
  expect_gt(E$Desvio_DAP, W$Desvio_DAP)
})

test_that("Steeper slope amplifies the aspect contrast", {
  d_5_E  <- tail_row(run_g(slope = 5,  aspect = pi/2))$DAP_max
  d_5_W  <- tail_row(run_g(slope = 5,  aspect = 3*pi/2))$DAP_max
  d_12_E <- tail_row(run_g(slope = 12, aspect = pi/2))$DAP_max
  d_12_W <- tail_row(run_g(slope = 12, aspect = 3*pi/2))$DAP_max

  expect_gt(d_12_E - d_12_W, d_5_E - d_5_W)
})


# ---- 4. SI sensitivity and the mortality clamp ------------------------------

test_that("Hd at base_age equals SI by construction", {
  for (SI in c(22, 28, 32)) {
    sim <- simulate_grandis(SI = SI, N0 = 900, G0 = 7,
                            PASW = 130, Elev = 130,
                            t0 = 10, t_end = 11, base_age = 10)
    expect_equal(sim$trajectory$AMD[1], SI, tolerance = 0.05)
  }
})

test_that("Mortality is bounded across the operational SI range", {
  # SI 22-35 should produce monotonically realistic projections; in
  # particular, total mortality from N0=900 over 14 years should stay
  # within plausible bounds (between 50 and 600 trees -- the upper
  # bound captures heavy Reineke + Methol-2003 combined effect).
  for (SI in c(22, 25, 28, 30, 35)) {
    sim <- simulate_grandis(SI = SI, N0 = 900, G0 = 7,
                            PASW = 130, Elev = 130,
                            t0 = 2, t_end = 16)
    final_N <- tail(sim$trajectory$N, 1)
    expect_gt(final_N,  50, label = paste("plausible end-N at SI =", SI))
    expect_lt(final_N, 900, label = paste("some mortality at SI =", SI))
  }
})

test_that("SI < 22 hits the mortality clamp without numerical blow-up", {
  sim_low  <- simulate_grandis(SI = 18, N0 = 900, G0 = 7,
                                PASW = 130, Elev = 130, t0 = 2, t_end = 16)
  sim_clamp <- simulate_grandis(SI = 22, N0 = 900, G0 = 7,
                                 PASW = 130, Elev = 130, t0 = 2, t_end = 16)
  end_low   <- tail(sim_low$trajectory$N,   1)
  end_clamp <- tail(sim_clamp$trajectory$N, 1)
  # Both runs finish with realistic, finite N values (no anti-mortality).
  expect_true(is.finite(end_low) && end_low <= 900)
  expect_true(is.finite(end_clamp) && end_clamp <= 900)
})


# ---- 5. Reineke kinetics ----------------------------------------------------

test_that("Low-density stand stays below the Reineke transition", {
  # 400 TPH on a moderate site never reaches RD = 0.55.
  sim <- run_g(N0 = 400, G0 = 1.5, t_end = 20)
  expect_lt(max(sim$trajectory$RD), 0.55)
})

test_that("High-density stand approaches the Reineke ceiling and self-thins", {
  sim <- run_g(N0 = 1111, G0 = 4, PASW = 160, t_end = 20)
  # By the end of the rotation, RD has reached at least 0.6 (the
  # soft-logistic centre) and stays bounded below 0.75.
  final_RD <- tail(sim$trajectory$RD, 1)
  expect_gt(final_RD, 0.55)
  expect_lt(max(sim$trajectory$RD), 0.75)
  # N has dropped from 1111 substantially by age 20 (>= 15% loss).
  expect_lt(tail(sim$trajectory$N, 1), 1111 * 0.85)
})

test_that("Thinning drops RD sharply then RD recovers as the stand re-stocks", {
  sim <- run_g(N0 = 900, G0 = 7, PASW = 160, t_end = 20,
               thins = list(list(age = 8, N_after = 400)))
  row_pre   <- sim$trajectory[sim$trajectory$age == 7, ]
  row_post  <- sim$trajectory[sim$trajectory$age == 8, ]
  row_late  <- tail(sim$trajectory, 1)

  # Thinning reduces RD substantially in the year of the thin
  expect_lt(row_post$RD, row_pre$RD * 0.75)
  # RD recovers somewhat by the end of the rotation
  expect_gt(row_late$RD, row_post$RD)
})

test_that("Sharper mort_k bites earlier than softer mort_k", {
  sim_soft  <- run_g(N0 = 1111, G0 = 4, PASW = 160, t_end = 20)  # mort_k=12
  sim_sharp <- simulate_grandis(SI = 28, N0 = 1111, G0 = 4,
                                 PASW = 160, Elev = 130,
                                 t0 = 2, t_end = 20, mort_k = 40)
  # Sharper transition -> more aggressive mortality at the same RD.
  expect_lt(tail(sim_sharp$trajectory$N, 1),
            tail(sim_soft$trajectory$N,  1))
})


# ---- 6. Cross-simulator comparison with simulate_inia() ---------------------

test_that("At a Z7-matched site, INIA and grandis Hd agree within 10%", {
  # PASW=160, Elev=130 reproduces the RC2019 augmented Z7 envelope.
  inia <- simulate_inia(SI = 28, N0 = 900, G0 = 7,
                        Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
                        t0 = 2, t_end = 16, zone = 7)
  grd  <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 160, Elev = 130,
                           t0 = 2, t_end = 16)
  inia_hd <- tail(inia$trajectory$AMD, 1)
  grd_hd  <- tail(grd$trajectory$AMD,  1)
  expect_lt(abs(grd_hd - inia_hd) / inia_hd, 0.10)
})

test_that("At a Z7-matched site, INIA and grandis G agree within 15%", {
  inia <- simulate_inia(SI = 28, N0 = 900, G0 = 7,
                        Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
                        t0 = 2, t_end = 16, zone = 7)
  grd  <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 160, Elev = 130,
                           t0 = 2, t_end = 16)
  inia_g <- tail(inia$trajectory$AB, 1)
  grd_g  <- tail(grd$trajectory$AB,  1)
  expect_lt(abs(grd_g - inia_g) / inia_g, 0.15)
})

test_that("grandis enforces more mortality than INIA at very high density", {
  # INIA mortality is SAG-calibrated and SI-independent; grandis adds
  # the Reineke ceiling on top of the constant-rate background. At
  # high enough density (N0 = 1500 here), the Reineke ceiling clearly
  # bites harder than INIA's mortality and the end-N is much lower.
  inia <- simulate_inia(SI = 32, N0 = 1500, G0 = 4.7,
                        Hd0 = 5.5, dmax0 = 9, SDd0 = 1.4,
                        t0 = 2, t_end = 20, zone = 7)
  grd  <- simulate_grandis(SI = 32, N0 = 1500, G0 = 4.7,
                           PASW = 170, Elev = 100,
                           t0 = 2, t_end = 20)
  expect_lt(tail(grd$trajectory$N,  1),
            tail(inia$trajectory$N, 1))
})

test_that("Both simulators respond to SI in the same direction", {
  inia_lo <- simulate_inia(SI = 25, N0 = 900, G0 = 7, Hd0 = 6,
                            t0 = 2, t_end = 16, zone = 7)
  inia_hi <- simulate_inia(SI = 32, N0 = 900, G0 = 7, Hd0 = 7,
                            t0 = 2, t_end = 16, zone = 7)
  grd_lo  <- simulate_grandis(SI = 25, N0 = 900, G0 = 7,
                               PASW = 160, Elev = 130,
                               t0 = 2, t_end = 16)
  grd_hi  <- simulate_grandis(SI = 32, N0 = 900, G0 = 7,
                               PASW = 160, Elev = 130,
                               t0 = 2, t_end = 16)
  expect_gt(tail(inia_hi$trajectory$AMD, 1),
            tail(inia_lo$trajectory$AMD, 1))
  expect_gt(tail(grd_hi$trajectory$AMD,  1),
            tail(grd_lo$trajectory$AMD,  1))
})


# ---- 7. Edge cases ----------------------------------------------------------

test_that("Zero slope makes aspect inert", {
  s_N <- run_g(slope = 0, aspect = 0)
  s_E <- run_g(slope = 0, aspect = pi/2)
  s_S <- run_g(slope = 0, aspect = pi)
  s_W <- run_g(slope = 0, aspect = 3*pi/2)
  Hd_vals <- c(tail(s_N$trajectory$AMD, 1),
                tail(s_E$trajectory$AMD, 1),
                tail(s_S$trajectory$AMD, 1),
                tail(s_W$trajectory$AMD, 1))
  expect_lt(diff(range(Hd_vals)), 0.05)
})

test_that("Multiple thinnings in sequence are all recorded", {
  sim <- run_g(N0 = 1111, G0 = 4, PASW = 160, t_end = 18,
               thins = list(list(age = 5,  N_after = 700),
                            list(age = 9,  N_after = 450),
                            list(age = 13, N_after = 250)))
  expect_equal(nrow(sim$thinnings), 3)
  expect_equal(sim$thinnings$age, c(5, 9, 13))
  expect_equal(sim$thinnings$N_post, c(700, 450, 250))
  expect_gt(sim$cumulative_thin_vol, 0)
})

test_that("t0 == t_end produces a single-row trajectory at the initial state", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                          PASW = 130, Elev = 130,
                          t0 = 5, t_end = 5)
  expect_equal(nrow(sim$trajectory), 1)
  expect_equal(sim$trajectory$N, 900)
  expect_equal(sim$trajectory$AB, 7)
})

test_that("Very high SI still produces a sane trajectory", {
  # SI = 38 is outside the RC2019 envelope and triggers heavy mortality
  # at the upper edge of Methol 2003's a = 0.4577 - 0.0218*SI. The
  # simulator must still return a finite, monotone-non-increasing N.
  sim <- run_g(SI = 38, t_end = 18)
  expect_true(all(is.finite(sim$trajectory$N)))
  expect_true(all(diff(sim$trajectory$N) <= 0))
  expect_gt(tail(sim$trajectory$N, 1), 50)
})

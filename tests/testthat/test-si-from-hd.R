# Tests for the inverse-Hd helpers: inia_si_from_hd / grandis_si_from_hd.

# ---- INIA ------------------------------------------------------------------

test_that("inia_si_from_hd recovers Hd at base_age (identity)", {
  # When age == base_age, the implied SI must equal the measured Hd.
  for (Hd in c(20, 25, 30, 35)) {
    for (zone in c(7, 8, 9)) {
      expect_equal(inia_si_from_hd(Hd, age = 10, zone = zone), Hd,
                   tolerance = 1e-8)
    }
  }
})

test_that("inia_si_from_hd round-trips with inia_hd_from_si", {
  # SI -> Hd at some age -> back to SI should be exact.
  for (SI in c(22, 28, 32)) {
    for (age in c(5, 8, 12, 16)) {
      Hd_at_age <- egrandis:::inia_hd_from_si(t = age, SI = SI, Z7 = 1)
      SI_back   <- inia_si_from_hd(Hd_at_age, age = age, zone = 7)
      expect_equal(SI_back, SI, tolerance = 1e-8,
                   info = sprintf("SI=%g, age=%g", SI, age))
    }
  }
})

test_that("inia_si_from_hd vectorises over Hd and age", {
  out <- inia_si_from_hd(Hd  = c(20, 24, 28),
                         age = c( 7,  8, 10), zone = 7)
  expect_length(out, 3)
  expect_true(all(out > 15) && all(out < 45))
})

test_that("inia_si_from_hd validates inputs", {
  expect_error(inia_si_from_hd(24, 8, zone = 6), "zone")
  expect_error(inia_si_from_hd(c(20, 24), age = 8), "same length")
  expect_error(inia_si_from_hd(-1, 8), "positive")
  expect_error(inia_si_from_hd(24, 0), "positive")
})


# ---- grandis ---------------------------------------------------------------

test_that("grandis_si_from_hd recovers Hd at base_age (identity)", {
  for (Hd in c(20, 25, 30)) {
    expect_equal(
      grandis_si_from_hd(Hd, age = 10, PASW = 130, slope = 0, aspect = 0),
      Hd, tolerance = 1e-8
    )
  }
})

test_that("grandis_si_from_hd round-trips with grandis_hd_from_si", {
  PASW <- 140; slope <- 6; aspect <- pi/4
  for (SI in c(22, 28, 32)) {
    for (age in c(5, 8, 12)) {
      Hd_at_age <- egrandis:::grandis_hd_from_si(
        t = age, SI = SI,
        PASW = PASW, slope = slope, aspect = aspect)
      SI_back <- grandis_si_from_hd(Hd_at_age, age = age,
                                     PASW = PASW, slope = slope, aspect = aspect)
      expect_equal(SI_back, SI, tolerance = 1e-8,
                   info = sprintf("SI=%g, age=%g", SI, age))
    }
  }
})

test_that("grandis_si_from_hd reflects PASW and aspect dependence", {
  # For age < base_age (the polymorphic projection going FORWARD from
  # the measurement to base_age), the curve at higher A (high PASW,
  # NE-facing) is shifted up at every age, so the SAME observed
  # (Hd, age) implies a HIGHER recovered SI -- not lower. This is the
  # round-trip behaviour: the recovered SI, fed back into the model
  # with the same site variables, reproduces the observation.
  base    <- grandis_si_from_hd(Hd = 24, age = 8,
                                 PASW = 100, slope = 0, aspect = 0)
  hi_pasw <- grandis_si_from_hd(Hd = 24, age = 8,
                                 PASW = 170, slope = 0, aspect = 0)
  expect_gt(hi_pasw, base)

  base_flat <- grandis_si_from_hd(Hd = 24, age = 8,
                                   PASW = 130, slope = 0, aspect = 0)
  ne_slope  <- grandis_si_from_hd(Hd = 24, age = 8,
                                   PASW = 130, slope = 9, aspect = pi/4)
  expect_gt(ne_slope, base_flat)

  # Past base_age the sign flips (the projection runs BACKWARD now), so
  # higher PASW implies LOWER recovered SI for the same (Hd, age).
  base_late    <- grandis_si_from_hd(Hd = 32, age = 14,
                                      PASW = 100, slope = 0, aspect = 0)
  hi_pasw_late <- grandis_si_from_hd(Hd = 32, age = 14,
                                      PASW = 170, slope = 0, aspect = 0)
  expect_lt(hi_pasw_late, base_late)
})

test_that("grandis_si_from_hd validates inputs", {
  expect_error(grandis_si_from_hd(c(20, 24), age = 8, PASW = 130),
               "same length")
  expect_error(grandis_si_from_hd(-1, 8, PASW = 130), "positive")
  expect_error(grandis_si_from_hd(24, 8, PASW = c(100, 130)), "scalars")
})

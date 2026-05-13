test_that("inia_tree_height saturates correctly", {
  # At d = dmax, height equals Hd
  expect_equal(inia_tree_height(49, Hd = 28, dmax = 49), 28, tolerance = 1e-6)
  expect_equal(inia_tree_height(35, Hd = 28, dmax = 35), 28, tolerance = 1e-6)

  # At d -> 0, height floor is 1.3 m (and the curve approaches 1.3 from above)
  expect_equal(inia_tree_height(0.1, Hd = 28, dmax = 49), 1.3, tolerance = 0.5)
  expect_gte(inia_tree_height(0.1, Hd = 28, dmax = 49), 1.3)

  # Monotonically increasing in d
  d <- seq(1, 50, length.out = 50)
  h <- inia_tree_height(d, Hd = 30, dmax = 55)
  expect_true(all(diff(h) >= 0))

  # Vectorised
  expect_length(inia_tree_height(c(10, 20, 30), Hd = 28, dmax = 49), 3)
})

test_that("inia_tree_agb scales as d^2.51 * h^0.24", {
  # Known sanity point
  expect_equal(inia_tree_agb(d = 35, h = 28), 665.1, tolerance = 1)

  # Doubling d should multiply AGB by 2^2.51 ~ 5.7
  ratio <- inia_tree_agb(50, 28) / inia_tree_agb(25, 28)
  expect_equal(ratio, 2^2.5094, tolerance = 0.01)

  # Doubling h should multiply AGB by 2^0.236 ~ 1.18
  ratio_h <- inia_tree_agb(30, 30) / inia_tree_agb(30, 15)
  expect_equal(ratio_h, 2^0.2362, tolerance = 0.01)
})

test_that("inia_stand_agb matches SAG 2021 reference points at age 10", {
  # Tolerances from handoff: +-5% for age 7+, +-15% for age 3-6.
  ref <- list(
    list(label = "Z7 SI=30 N=550 age 10",
         args = list(N = 364, Dq = 38.3, dmax = 49.2, SDd = 6.1, Hd = 30.0),
         target = 306.5),
    list(label = "Z7 SI=25 N=550 age 10",
         args = list(N = 364, Dq = 38.3, dmax = 49.2, SDd = 6.1, Hd = 25.0),
         target = 286.1),
    list(label = "Z7 SI=30 N=1111 age 10",
         args = list(N = 737, Dq = 27.7, dmax = 42.9, SDd = 6.1, Hd = 30.0),
         target = 290.0)
  )
  for (r in ref) {
    got <- do.call(inia_stand_agb, r$args)$agb
    rel <- abs(got - r$target) / r$target
    expect_lt(rel, 0.05, label = r$label)
  }
})

test_that("inia_stand_agb returns 0 for empty / degenerate stands", {
  out <- inia_stand_agb(N = 0, Dq = 0, dmax = 0.1, SDd = 0.01, Hd = 1.4)
  expect_equal(out$agb, 0)
  expect_equal(out$carbon, 0)
  expect_equal(out$co2eq, 0)
})

test_that("inia_add_biomass appends exactly 3 columns and preserves others", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  cols_before <- names(sim$trajectory)
  sim2 <- inia_add_biomass(sim)
  cols_after <- names(sim2$trajectory)
  new_cols <- setdiff(cols_after, cols_before)
  expect_setequal(new_cols, c("Biomasa", "Carbon", "CO2eq"))

  # Original columns unchanged
  for (col in cols_before) {
    expect_equal(sim2$trajectory[[col]], sim$trajectory[[col]],
                 label = paste("column", col))
  }

  # Non-trajectory parts of the sim object are unchanged.
  for (nm in setdiff(names(sim), "trajectory")) {
    expect_equal(sim2[[nm]], sim[[nm]])
  }
})

test_that("CO2eq is consistently 1.797 * Biomasa across the trajectory", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  sim <- inia_add_biomass(sim)
  ratios <- sim$trajectory$CO2eq / sim$trajectory$Biomasa
  # The ratios should be tight around 1.797; rounding to 1 decimal place
  # introduces small noise so we allow a 0.005 band.
  expect_true(all(abs(ratios - 1.797) < 0.005))
})

test_that("Carbon is consistently 0.49 * Biomasa across the trajectory", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  sim <- inia_add_biomass(sim)
  # At very low ages B is small enough that 1-decimal rounding distorts
  # the ratio (e.g. B=3.4, C=1.7, ratio=0.500 instead of 0.49). Use a
  # 0.011 tolerance to cover the worst-case rounding swing.
  ratios <- sim$trajectory$Carbon / sim$trajectory$Biomasa
  expect_true(all(abs(ratios - 0.49) < 0.011))
})

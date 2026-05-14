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

test_that("tree-level components match Winck 2015 Table 3 equations", {
  # Direct evaluation of the Zone 1 NE Argentina equations.
  d <- 22.3; h <- 26.7

  bt_ref <- 1.01 * exp(-3.36 + 2.12 * log(d) +  0.65 * log(h))
  bf_ref <- 1.01 * exp(-4.51 + 1.83 * log(d) +  1.22 * log(h))
  br_ref <- 1.13 * exp(-2.68 + 3.73 * log(d) + -1.77 * log(h))

  expect_equal(inia_tree_agb(d, h),      bt_ref, tolerance = 1e-6)
  expect_equal(inia_tree_stem(d, h),     bf_ref, tolerance = 1e-6)
  expect_equal(inia_tree_branches(d, h), br_ref, tolerance = 1e-6)
})

test_that("tree-level helpers are vectorised and monotone in d", {
  d <- c(15, 25, 35)
  h <- c(18, 26, 30)
  expect_length(inia_tree_stem(d, h),     3)
  expect_length(inia_tree_branches(d, h), 3)
  expect_length(inia_tree_agb(d, h),      3)
  expect_true(all(inia_tree_agb(d, h) > 0))

  d_grid <- seq(12, 36, length.out = 50)
  h_grid <- rep(26.7, 50)
  expect_true(all(diff(inia_tree_agb(d_grid, h_grid)) > 0))
})

test_that("Winck per-tree AGB lands within physical wood-density bounds", {
  # Wood + bark dry density for E. grandis is bounded by ~600 kg/m^3
  # in the published literature (Resquin 2019 Uruguay, Vital 1984
  # Brazil). The fitted regression residual lifts the implied AGB
  # density slightly above pure stem density at small trees (more
  # crown fraction), so we cap the practical envelope at 650 kg/m^3.
  V_fang <- function(d, h) 4e-5 * d^2.09 * h^0.862
  test_pts <- list(c(15, 18), c(22, 25), c(28, 30), c(35, 40))
  for (p in test_pts) {
    d <- p[1]; h <- p[2]
    rho <- inia_tree_agb(d, h) / V_fang(d, h)
    expect_lt(rho, 650, label = sprintf("rho_implied at d=%g,h=%g", d, h))
    expect_gt(rho, 350, label = sprintf("rho_implied at d=%g,h=%g", d, h))
  }
})

test_that("Winck stem fraction tracks the source paper average (~85%)", {
  # Winck Table 2 zone 1: stem 85.8%, branches 11.5%, leaves 2.6%.
  # Components and AGB were fit independently, so allow a wide window.
  d_grid <- c(20, 25, 30, 35)
  h_grid <- c(22, 28, 32, 38)
  for (i in seq_along(d_grid)) {
    frac <- inia_tree_stem(d_grid[i], h_grid[i]) /
            inia_tree_agb(d_grid[i],  h_grid[i])
    expect_gt(frac, 0.70)
    expect_lt(frac, 0.95)
  }
})

test_that("inia_stand_agb returns stem and branches alongside agb", {
  out <- inia_stand_agb(N = 550, Dq = 24.0, dmax = 36.0,
                        SDd = 4.5, Hd = 27.0)
  expect_true(all(c("agb", "stem", "branches", "carbon", "co2eq") %in% names(out)))
  expect_gt(out$stem, 0)
  expect_gt(out$branches, 0)
  # Stem dominates the residual; branches are non-negligible but smaller.
  expect_gt(out$stem, out$branches)
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
  # 1-decimal rounding of small Biomasa values (age 1, Biomasa ~ 3 t/ha)
  # makes the ratio noisy, so skip ages where the rounding artefact
  # dominates and check the rest tightly.
  keep <- sim$trajectory$Biomasa >= 5
  ratios <- sim$trajectory$CO2eq[keep] / sim$trajectory$Biomasa[keep]
  expect_true(all(abs(ratios - 1.797) < 0.005))
})

test_that("Carbon is consistently 0.49 * Biomasa across the trajectory", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  sim <- inia_add_biomass(sim)
  # 1-decimal rounding of small Biomasa values distorts the ratio;
  # skip the rounding-dominated rows for the same reason as the CO2eq
  # test above.
  keep <- sim$trajectory$Biomasa >= 5
  ratios <- sim$trajectory$Carbon[keep] / sim$trajectory$Biomasa[keep]
  expect_true(all(abs(ratios - 0.49) < 0.005))
})

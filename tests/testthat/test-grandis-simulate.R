# Phase 2 tests for R/grandis_simulate.R.
# Sensitivity / Reineke-trajectory tests live alongside in Phase 3.

# ---- shape and provenance ---------------------------------------------------

test_that("simulate_grandis returns the simulate_inia trajectory schema", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 130, Elev = 130,
                           t0 = 2, t_end = 12)
  expect_named(sim, c("trajectory", "thinnings", "cumulative_thin_vol",
                      "final_standing", "total_yield", "parameters",
                      "provenance"))
  expect_setequal(
    colnames(sim$trajectory),
    c("age","AMD","N","AB","DAP_medio","DAP_max","Desvio_DAP",
      "Vol_Total","MAI","SDI","RD","ICA")
  )
  expect_equal(nrow(sim$trajectory), 11)
  expect_null(sim$thinnings)
})

test_that("provenance object has the expected structure", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 130, Elev = 130, t0 = 2, t_end = 5)
  p <- sim$provenance
  expect_equal(p$model, "simulate_grandis")
  expect_true(is.character(p$caveats) && length(p$caveats) >= 4)
  expect_true(all(c("Hd", "G", "dmax", "SDd", "N_bg", "N_self_thin", "V")
                   %in% names(p$submodels)))
  expect_true(all(c("SI", "PASW", "Elev", "slope", "aspect")
                   %in% names(p$inputs)))
})

# ---- accumulated invariants -------------------------------------------------

test_that("MAI = Vol_Total / age for every trajectory row", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 130, Elev = 130,
                           t0 = 2, t_end = 16)
  # The trajectory rounds independently rather than computing MAI from
  # the rounded Vol_Total, so half-to-even rounding can produce 0.1
  # ulp discrepancies. Compare on the un-rounded ratio.
  mai_ref <- sim$trajectory$Vol_Total / sim$trajectory$age
  expect_lt(max(abs(sim$trajectory$MAI - mai_ref)), 0.1)
})

test_that("ICA equals the row-to-row difference in Vol_Total", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 130, Elev = 130,
                           t0 = 2, t_end = 16)
  diffs <- c(sim$trajectory$MAI[1],
             round(diff(sim$trajectory$Vol_Total), 1))
  expect_equal(sim$trajectory$ICA, diffs)
})

test_that("N is non-increasing across the trajectory", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 130, Elev = 130,
                           t0 = 2, t_end = 20)
  expect_true(all(diff(sim$trajectory$N) <= 0))
})

# ---- input validation -------------------------------------------------------

test_that("simulate_grandis rejects negative PASW, Elev, slope", {
  expect_error(simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                                 PASW = -1, Elev = 130), "PASW")
  expect_error(simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                                 PASW = 130, Elev = -1), "Elev")
  expect_error(simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                                 PASW = 130, Elev = 130, slope = -5), "slope")
})

# ---- thinning ---------------------------------------------------------------

test_that("a single thinning event is recorded and modifies the trajectory", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 130, Elev = 130,
                           t0 = 2, t_end = 14,
                           thins = list(list(age = 7, N_after = 500)))
  expect_equal(nrow(sim$thinnings), 1)
  expect_equal(sim$thinnings$age, 7)
  expect_equal(sim$thinnings$N_post, 500)
  expect_gt(sim$thinnings$V_removed, 0)
  expect_gt(sim$cumulative_thin_vol, 0)
  # N drops sharply at the thin year
  row_pre  <- sim$trajectory[sim$trajectory$age == 6, ]
  row_post <- sim$trajectory[sim$trajectory$age == 7, ]
  expect_gt(row_pre$N, row_post$N)
})

# ---- site-variable sensitivity ----------------------------------------------

test_that("higher PASW produces taller, more productive stands", {
  sim_lo <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                              PASW = 80, Elev = 130, t0 = 2, t_end = 16)
  sim_hi <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                              PASW = 160, Elev = 130, t0 = 2, t_end = 16)
  end_lo <- tail(sim_lo$trajectory, 1)
  end_hi <- tail(sim_hi$trajectory, 1)
  expect_gt(end_hi$AMD,       end_lo$AMD)
  expect_gt(end_hi$AB,        end_lo$AB)
  expect_gt(end_hi$Vol_Total, end_lo$Vol_Total)
})

test_that("Reineke self-thinning bounds high-density rotations", {
  # Two stands starting at 400 vs 1100 TPH on the same site should
  # converge in RD over time, with the high-density stand experiencing
  # substantially more mortality.
  s_lo <- simulate_grandis(SI = 28, N0 =  400, G0 = 1.5,
                            PASW = 140, Elev = 130, t0 = 2, t_end = 18)
  s_hi <- simulate_grandis(SI = 28, N0 = 1100, G0 = 4.0,
                            PASW = 140, Elev = 130, t0 = 2, t_end = 18)
  end_lo <- tail(s_lo$trajectory, 1)
  end_hi <- tail(s_hi$trajectory, 1)
  # High-density stand kept more trees but is closer to its RD ceiling.
  expect_gt(end_hi$N, end_lo$N)
  expect_gt(end_hi$RD, end_lo$RD)
  # Both stands settled below RD = 0.75 (the soft logistic flattens here).
  expect_lt(end_hi$RD, 0.75)
})

# ---- downstream-tool compatibility ------------------------------------------

test_that("downstream tools (merch_vol, add_biomass, dmd_plot) accept the result", {
  sim <- simulate_grandis(SI = 28, N0 = 900, G0 = 7,
                           PASW = 130, Elev = 130,
                           t0 = 2, t_end = 16)
  expect_no_error(inia_add_biomass(sim))
  expect_no_error(inia_merch_vol(sim, age = 12))
  expect_no_error(inia_get_distribution(sim, age = 10))

  skip_if_not_installed("ggplot2")
  expect_no_error(inia_dmd_plot(sim))
})

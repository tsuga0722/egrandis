test_that("inia_dmd_plot returns a ggplot object on an unthinned trajectory", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_inia(SI = 28, N0 = 900, G0 = 7,
                       t0 = 3, t_end = 15, zone = 7)
  p <- inia_dmd_plot(sim)
  expect_s3_class(p, "ggplot")
  # Reference lines + trajectory path + trajectory points + age labels = 4.
  # Thinnings would add a 5th. Allow >= 4 to be forward-compatible.
  expect_gte(length(p$layers), 4)
})

test_that("inia_dmd_plot adds a thinning-arrow layer when thins are present", {
  skip_if_not_installed("ggplot2")
  sim_un  <- simulate_inia(SI = 28, N0 = 900, G0 = 7,
                           t0 = 3, t_end = 15, zone = 7)
  sim_thn <- simulate_inia(SI = 28, N0 = 900, G0 = 7,
                           t0 = 3, t_end = 15, zone = 7,
                           thins = list(list(age = 7, N_after = 500)))
  expect_equal(length(inia_dmd_plot(sim_thn)$layers),
               length(inia_dmd_plot(sim_un)$layers) + 1)
})

test_that("inia_dmd_plot respects custom axis ranges", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_inia(SI = 28, N0 = 900, G0 = 7,
                       t0 = 3, t_end = 15, zone = 7)
  p <- inia_dmd_plot(sim, dq_range = c(8, 50), n_range = c(200, 2000))
  # coord_cartesian stores limits on the coord object
  expect_equal(p$coordinates$limits$x, c(8, 50))
  expect_equal(p$coordinates$limits$y, c(200, 2000))
})

test_that("Reineke reference line satisfies N(Dq) = SDImax * (25/Dq)^1.605", {
  # Closed-form check that the function's parameters are wired correctly.
  p   <- .inia_params$sdi
  Dq  <- 25
  N_max_at_25 <- p$SDImax * (p$Dq_ref / Dq)^p$beta
  expect_equal(N_max_at_25, p$SDImax)  # at Dq = Dq_ref, max line = SDImax

  Dq  <- 10
  N_expected <- 1250 * (25 / 10)^1.605
  expect_equal(p$SDImax * (p$Dq_ref / Dq)^p$beta, N_expected)
})

test_that("RD column equals SDI / SDImax in the trajectory", {
  sim <- simulate_inia(SI = 28, N0 = 900, G0 = 7,
                       t0 = 3, t_end = 15, zone = 7)
  expected <- round(sim$trajectory$SDI / .inia_params$sdi$SDImax, 3)
  # Allow a tiny tolerance because SDI is itself rounded for display.
  expect_lt(max(abs(sim$trajectory$RD - expected)), 0.002)
  expect_true(all(sim$trajectory$RD >= 0 & sim$trajectory$RD <= 1.5))
})

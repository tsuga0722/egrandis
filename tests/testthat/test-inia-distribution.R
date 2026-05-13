test_that("Weibull recovery returns sensible parameters", {
  p <- egrandis:::inia_weibull_params(Dq = 38.3, dmax = 49.2, SDd = 6.1)
  expect_equal(p$a, 49.2)                # location is dmax
  expect_gt(p$b, 0)
  expect_gt(p$c, 1)
})

test_that("Weibull CDF is monotonic on [0, dmax]", {
  pp <- egrandis:::inia_weibull_params(Dq = 30, dmax = 45, SDd = 5)
  x <- seq(0.1, 44.9, length.out = 200)
  cdf <- vapply(x, egrandis:::inia_weibull_cdf, numeric(1),
                a = pp$a, b = pp$b, c = pp$c)
  expect_true(all(diff(cdf) >= 0))
  expect_equal(egrandis:::inia_weibull_cdf(45, pp$a, pp$b, pp$c), 1)
  expect_equal(egrandis:::inia_weibull_cdf(0,  pp$a, pp$b, pp$c), 0)
})

test_that("inia_diam_dist sums to (approx) N and yields BA close to G", {
  N  <- 364
  Dq <- 38.3
  G_expected <- N * pi * (Dq / 200)^2  # ~41.9 m2/ha
  dd <- inia_diam_dist(N = N, Dq = Dq, dmax = 49.2, SDd = 6.1)
  expect_equal(sum(dd$freq), N, tolerance = 0.05)
  expect_equal(sum(dd$ba_ha), G_expected, tolerance = 0.05)
})

test_that("inia_get_distribution pulls the right row from a simulation", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  dd10 <- inia_get_distribution(sim, age = 10)
  expect_s3_class(dd10, "data.frame")
  expect_true(all(c("class_mid", "freq", "ba_ha") %in% names(dd10)))
  expect_error(inia_get_distribution(sim, age = 99), "not found")
})

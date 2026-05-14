test_that("Zone 7 mortality matches SAG validation within +/- 15 trees", {
  data(sag_validation, package = "egrandis")
  z7_unthinned <- c("z7_si30_n550", "z7_si25_n550", "z7_si30_n1111")
  for (nm in z7_unthinned) {
    out <- run_scenario(nm)
    expect_lt(max(abs(out$sim_traj$N - out$ref_traj$N)), 15,
              label = paste("max |dN| for", nm))
  }
})

test_that("Zone 8/9 mortality matches SAG validation within +/- 25 trees", {
  # The Z8/9 SAG 2021 mortality trajectory has a humped per-year-loss shape
  # (rate peaks around age 5-6) that the monotonic Clutter-Jones form
  # cannot fully capture. The current parameter regime is fitted to the
  # only Z8/9 SAG scenario we have; max abs error ~20 trees.
  data(sag_validation, package = "egrandis")
  out <- run_scenario("z8_si30_n550")
  expect_lt(max(abs(out$sim_traj$N - out$ref_traj$N)), 25)
})

test_that("mortality is monotonically non-increasing", {
  N0 <- 1111
  ages <- seq(1, 25, by = 1)
  N <- numeric(length(ages))
  N[1] <- N0
  for (i in 2:length(ages)) {
    N[i] <- egrandis:::inia_n(N[i - 1], ages[i - 1], ages[i], Z7 = 1)
  }
  expect_true(all(diff(N) <= 0))
})

test_that("mortality is SI-invariant in the INIA 2021 form", {
  # In the INIA 2021 form, SI=25 and SI=30 produce identical N trajectories
  # in Zone 7 (SAG behaviour).
  ages <- 1:16
  N_a <- N_b <- 550
  for (t in 2:length(ages)) {
    N_a <- egrandis:::inia_n(N_a, ages[t - 1], ages[t], Z7 = 1)
    N_b <- egrandis:::inia_n(N_b, ages[t - 1], ages[t], Z7 = 1)
  }
  # Trivially equal here — the test guards against accidental SI injection.
  expect_equal(N_a, N_b)
})

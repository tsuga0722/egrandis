test_that("dominant height closes at the base age (SI invariance)", {
  for (SI in c(20, 25, 30, 35)) {
    expect_equal(egrandis:::inia_hd_from_si(10, SI = SI, Z7 = 1), SI,
                 tolerance = 1e-9, info = paste("Z7, SI =", SI))
    expect_equal(egrandis:::inia_hd_from_si(10, SI = SI, Z7 = 0), SI,
                 tolerance = 1e-9, info = paste("Z8/9, SI =", SI))
  }
})

test_that("Hd is monotonically increasing with age", {
  ages <- 1:30
  hd <- vapply(ages, egrandis:::inia_hd_from_si, numeric(1), SI = 30, Z7 = 1)
  expect_true(all(diff(hd) > 0))
})

test_that("higher Z7 asymptote means higher Hd at any age (post-base-age)", {
  ages <- 11:25
  hd_z7  <- vapply(ages, egrandis:::inia_hd_from_si, numeric(1), SI = 30, Z7 = 1)
  hd_z89 <- vapply(ages, egrandis:::inia_hd_from_si, numeric(1), SI = 30, Z7 = 0)
  # Same SI at base age 10, but Zone 7 has the higher asymptote, so post-10
  # Zone 7 should track higher than Zone 8/9.
  expect_true(all(hd_z7 >= hd_z89))
})

test_that("Hd matches SAG validation trajectories for all unthinned scenarios", {
  data(sag_validation, package = "egrandis")
  for (nm in setdiff(names(sag_validation), "z7_si30_n550_thinned")) {
    out <- run_scenario(nm)
    expect_lt(max(abs(out$sim_traj$AMD - out$ref_traj$AMD)), 0.2,
              label = paste("max |dHd| for", nm))
  }
})

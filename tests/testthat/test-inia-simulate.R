test_that("simulate_inia returns the expected shape", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  expect_named(sim, c("trajectory", "thinnings", "cumulative_thin_vol",
                      "final_standing", "total_yield", "parameters"))
  expect_s3_class(sim$trajectory, "data.frame")
  expect_equal(nrow(sim$trajectory), 16)
  expect_true(all(c("age", "AMD", "N", "AB", "DAP_medio", "DAP_max",
                    "Desvio_DAP", "Vol_Total", "MAI", "SDI", "ICA")
                  %in% names(sim$trajectory)))
  expect_null(sim$thinnings)
  expect_equal(sim$cumulative_thin_vol, 0)
})

test_that("simulate_inia rejects unknown zones", {
  expect_error(
    simulate_inia(SI = 30, N0 = 550, G0 = 1.7,
                  Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
                  t0 = 1, t_end = 5, zone = 0),
    "must be 7, 8, or 9"
  )
})

test_that("thinnings table records each thin event in order", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.16, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 12, zone = 7,
    thins = list(list(age = 3, N_after = 412),
                 list(age = 7, N_after = 197))
  )
  expect_s3_class(sim$thinnings, "data.frame")
  expect_equal(nrow(sim$thinnings), 2)
  expect_equal(sim$thinnings$age, c(3, 7))
  expect_equal(sim$thinnings$N_post, c(412, 197))
  expect_true(sim$cumulative_thin_vol > 0)
  expect_equal(
    sim$total_yield,
    round(sim$final_standing + sim$cumulative_thin_vol, 1),
    tolerance = 0.2
  )
})

test_that("integration: unthinned z7_si30_n550 reproduces SAG trajectory", {
  data(sag_validation, package = "egrandis")
  out <- run_scenario("z7_si30_n550")
  sim <- out$sim_traj; ref <- out$ref_traj
  # Documented per-submodel tolerances (see README "Model Details").
  expect_lt(max(abs(sim$AMD - ref$AMD)),  0.2)
  expect_lt(max(abs(sim$AB  - ref$AB )),  0.5)
  expect_lt(max(abs(sim$N   - ref$N  )),  15)
  expect_lt(max(abs(sim$Desvio_DAP - ref$Desvio_DAP)), 0.1)
  expect_lt(max(abs(sim$DAP_max - ref$DAP_max)), 1.5)
})

test_that("integration: high density (1111 TPH) trajectory tracks SAG", {
  # The dmax submodel was fit only at Zone 7 / 550 TPH.
  # We do not assert dmax here because SAG's high-density trajectory diverges.
  data(sag_validation, package = "egrandis")
  out <- run_scenario("z7_si30_n1111")
  sim <- out$sim_traj; ref <- out$ref_traj
  expect_lt(max(abs(sim$AMD - ref$AMD)), 0.2)
  expect_lt(max(abs(sim$AB  - ref$AB )), 0.5)
  expect_lt(max(abs(sim$N   - ref$N  )), 15)
})

test_that("integration: Zone 8 trajectory tracks SAG", {
  data(sag_validation, package = "egrandis")
  out <- run_scenario("z8_si30_n550")
  sim <- out$sim_traj; ref <- out$ref_traj
  expect_lt(max(abs(sim$AMD - ref$AMD)),  0.2)
  expect_lt(max(abs(sim$AB  - ref$AB )),  0.6)
  expect_lt(max(abs(sim$N   - ref$N  )), 25)
  expect_lt(max(abs(sim$DAP_max - ref$DAP_max)), 1.0)
})

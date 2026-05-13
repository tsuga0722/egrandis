test_that("basal area matches SAG validation for unthinned scenarios", {
  data(sag_validation, package = "egrandis")
  for (nm in setdiff(names(sag_validation), "z7_si30_n550_thinned")) {
    out <- run_scenario(nm)
    expect_lt(max(abs(out$sim_traj$AB - out$ref_traj$AB)), 0.5,
              label = paste("max |dG| for", nm))
  }
})

test_that("basal area asymptotes to the Z7 / Z8-9 limits", {
  # Project a large stand forward many years; G should approach
  # exp(a0 + a1*Z7).
  G_z7  <- egrandis:::inia_g(40, 10, 200, Z7 = 1)
  G_z89 <- egrandis:::inia_g(40, 10, 200, Z7 = 0)
  expect_equal(G_z7,  exp(4.027), tolerance = 0.5)
  expect_equal(G_z89, exp(3.753), tolerance = 0.5)
})

test_that("BA projection without thinning is a no-op of itself (idempotent reprojection)", {
  # Projecting t1 -> t2 then t2 -> t2 should equal t1 -> t2.
  G_a <- egrandis:::inia_g(20, 5, 10, Z7 = 1)
  G_b <- egrandis:::inia_g(G_a, 10, 10, Z7 = 1)
  expect_equal(G_a, G_b, tolerance = 1e-9)
})

test_that("thinning modifier reduces the C exponent (slows convergence)", {
  G1 <- 25
  # Without thinning history
  G_unthinned <- egrandis:::inia_g(G1, 7, 12, Z7 = 1)
  # With recent thinning: Na/Nb = 0.5 at age 7
  G_thinned <- egrandis:::inia_g(G1, 7, 12, Z7 = 1, Na = 250, Nb = 500, tt = 7)
  # Thinned stand should re-occupy slower -> lower G at the target age
  expect_lt(G_thinned, G_unthinned)
})

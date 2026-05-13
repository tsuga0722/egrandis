test_that("taper at breast height returns approximately DBH", {
  # Realistic D-H combinations within the Hirigoyen 2021 calibration
  # envelope (D 11-37 cm, H 17-41 m). Outside this range slenderness
  # ratios become unrealistic and the Fang fit drifts >5%.
  cases <- list(
    c(15, 18), c(15, 25),
    c(22, 22), c(22, 28), c(22, 35),
    c(30, 28), c(30, 35),
    c(40, 35)
  )
  for (cs in cases) {
    D <- cs[1]; H <- cs[2]
    d_bh <- inia_taper(1.3, D = D, H = H)
    expect_lt(abs(d_bh - D) / D, 0.05,
              label = sprintf("d at 1.3m for D=%d, H=%d", D, H))
  }
})

test_that("taper is 0 at the tip and below the stump", {
  expect_equal(inia_taper(28,   D = 25, H = 28), 0)
  expect_equal(inia_taper(0.05, D = 25, H = 28), 0)  # below default stump
  expect_equal(inia_taper(50,   D = 25, H = 28), 0)  # above tip
})

test_that("taper is monotonically non-increasing from stump to tip", {
  h <- seq(0.10, 28, length.out = 300)
  d <- inia_taper(h, D = 25, H = 28)
  expect_true(all(diff(d) <= 1e-6))
})

test_that("integrated tree volume matches the closed-form Fang total", {
  # Cone tail above p2 loses < 1% of the formula's segment-3 volume,
  # which is itself < 1% of the tree.
  for (D in c(15, 25, 35)) {
    for (H in c(15, 25, 35)) {
      Vc <- inia_tree_total_vol(D, H)
      Vi <- inia_tree_vol(D, H, h_lower = 0.10, h_upper = H)
      expect_lt(abs(Vi - Vc) / Vc, 0.01,
                label = sprintf("V_closed vs V_int for D=%d H=%d", D, H))
    }
  }
})

test_that("merchantable volume decreases monotonically as small-end limit grows", {
  D <- 30; H <- 32
  V_total <- inia_tree_total_vol(D, H)
  V8  <- inia_tree_vol(D, H, d_top = 8)
  V14 <- inia_tree_vol(D, H, d_top = 14)
  V25 <- inia_tree_vol(D, H, d_top = 25)
  expect_lt(V25, V14)
  expect_lt(V14, V8)
  expect_lt(V8,  V_total)
})

test_that("inia_height_at_d inverts the taper", {
  D <- 30; H <- 32
  for (d_top in c(8, 14, 20, 25)) {
    h <- inia_height_at_d(d_top, D, H)
    d_back <- inia_taper(h, D, H)
    expect_lt(abs(d_back - d_top), 0.05)
  }
})

test_that("inia_height_at_d handles edge cases", {
  expect_equal(inia_height_at_d(d_top = 0,  D = 25, H = 28), 28)  # never reaches 0
  expect_equal(inia_height_at_d(d_top = 99, D = 25, H = 28), 0.10) # always above 99cm? No: butt < 99
})

test_that("inia_height_class produces sensible per-class heights", {
  h <- inia_height_class(D = c(10, 30, 50), Dq = 30, Hd = 28)
  expect_length(h, 3)
  expect_true(all(h > 1.3))
  expect_true(h[1] < h[2] && h[2] <= h[3])
  # Dq tree should hit ~Hd
  expect_equal(inia_height_class(30, Dq = 30, Hd = 28), 28, tolerance = 0.01)
})

test_that("inia_merch_vol returns shape and sane totals", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  mv <- inia_merch_vol(sim, age = 16)
  expect_s3_class(mv, "data.frame")
  expect_true(all(c("class_mid", "freq", "H_class", "vol_total",
                    "vol_veneer", "vol_solid", "vol_pulp",
                    "vol_top_waste") %in% names(mv)))

  totals <- attr(mv, "totals")
  # Sum of per-product volumes plus top waste should equal total volume.
  parts_sum <- totals["vol_veneer"] + totals["vol_solid"] +
               totals["vol_pulp"]   + totals["top_waste"]
  expect_equal(unname(parts_sum), unname(totals["total"]),
               tolerance = 0.01)

  # All product volumes are non-negative and finite.
  expect_true(all(is.finite(unname(totals))))
  expect_true(all(unname(totals) >= 0))
})

test_that("inia_merch_vol respects custom product specs", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 16, zone = 7
  )
  mv <- inia_merch_vol(
    sim, age = 16,
    products = list(saw = list(d_min = 18), pulp = list(d_min = 8))
  )
  totals <- attr(mv, "totals")
  expect_true("vol_saw"  %in% names(totals))
  expect_true("vol_pulp" %in% names(totals))
  expect_false("vol_veneer" %in% names(totals))
})

test_that("inia_merch_vol errors on missing age", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 12, zone = 7
  )
  expect_error(inia_merch_vol(sim, age = 99), "not found")
})

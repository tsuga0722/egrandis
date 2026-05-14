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
                    "vol_large_sawlog", "vol_small_sawlog", "vol_pulp",
                    "vol_top_waste") %in% names(mv)))

  totals <- attr(mv, "totals")
  # Sum of per-product volumes plus top waste should equal total volume.
  parts_sum <- totals["vol_large_sawlog"] + totals["vol_small_sawlog"] +
               totals["vol_pulp"]         + totals["top_waste"]
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
  # Default product names should not leak in when a custom list is passed.
  expect_false("vol_large_sawlog" %in% names(totals))
})

test_that("inia_merch_vol errors on missing age", {
  sim <- simulate_inia(
    SI = 30, N0 = 550, G0 = 1.7,
    Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
    t0 = 1, t_end = 12, zone = 7
  )
  expect_error(inia_merch_vol(sim, age = 99), "not found")
})

# ---- Min log length ---------------------------------------------------------

# Reusable fixture for the new-feature tests.
make_mid_rotation_sim <- function() {
  simulate_inia(
    SI = 28, N0 = 900, G0 = 7.0,
    Hd0 = 7.0, dmax0 = 13.0, SDd0 = 1.8,
    t0 = 2, t_end = 16, zone = 7
  )
}

test_that("l_min defaults to no constraint and matches prior behavior", {
  sim <- make_mid_rotation_sim()
  a <- attr(inia_merch_vol(sim, age = 16), "totals")
  b <- attr(inia_merch_vol(sim, age = 16,
              products = list(
                large_sawlog = list(d_min = 25),
                small_sawlog = list(d_min = 14),
                pulp         = list(d_min = 8)
              )), "totals")
  expect_equal(unname(a["vol_large_sawlog"]), unname(b["vol_large_sawlog"]),
               tolerance = 1e-6)
  expect_equal(unname(a["total"]), unname(b["total"]), tolerance = 1e-6)
})

test_that("a binding l_min shifts volume from a higher grade to a lower one", {
  sim <- make_mid_rotation_sim()
  unconstr <- attr(inia_merch_vol(sim, age = 16), "totals")
  constrained <- attr(inia_merch_vol(sim, age = 16,
    products = list(
      large_sawlog = list(d_min = 25, l_min = 4),
      small_sawlog = list(d_min = 14, l_min = 2.4),
      pulp         = list(d_min = 8)
    )), "totals")
  # Large sawlog volume must drop, small sawlog must rise, and the
  # stand-level total bole volume is conserved.
  expect_lt(constrained[["vol_large_sawlog"]], unconstr[["vol_large_sawlog"]])
  expect_gt(constrained[["vol_small_sawlog"]], unconstr[["vol_small_sawlog"]])
  expect_equal(constrained[["total"]], unconstr[["total"]], tolerance = 1e-6)
})

test_that("an unreachable l_min skips the grade entirely", {
  sim <- make_mid_rotation_sim()
  mv <- inia_merch_vol(sim, age = 16,
    products = list(
      large_sawlog = list(d_min = 25, l_min = 100),  # nothing is 100 m long
      small_sawlog = list(d_min = 14),
      pulp         = list(d_min = 8)
    ))
  totals <- attr(mv, "totals")
  expect_equal(totals[["vol_large_sawlog"]], 0)
  # All bole volume now falls to small_sawlog + pulp + top_waste; total
  # bole conserved.
  expect_equal(totals[["total"]],
               attr(inia_merch_vol(sim, age = 16), "totals")[["total"]],
               tolerance = 1e-6)
})

# ---- Pruned-height split ----------------------------------------------------

test_that("pruned_height NULL behaves identically to the unsplit call", {
  sim <- make_mid_rotation_sim()
  mv_null <- inia_merch_vol(sim, age = 16, pruned_height = NULL)
  mv_base <- inia_merch_vol(sim, age = 16)
  expect_equal(attr(mv_null, "totals"), attr(mv_base, "totals"))
  expect_false(any(grepl("_pruned$",   names(mv_null))))
  expect_false(any(grepl("_unpruned$", names(mv_null))))
})

test_that("pruned + unpruned columns sum to the parent product volume", {
  sim <- make_mid_rotation_sim()
  mv <- inia_merch_vol(sim, age = 16, pruned_height = 6)
  for (nm in c("vol_large_sawlog", "vol_small_sawlog", "vol_pulp")) {
    parent <- attr(mv, "totals")[[nm]]
    split  <- attr(mv, "totals")[[paste0(nm, "_pruned")]] +
              attr(mv, "totals")[[paste0(nm, "_unpruned")]]
    expect_equal(parent, split, tolerance = 1e-6,
                 info = paste("product:", nm))
  }
})

test_that("pruned_height above the tallest tree pushes all volume to pruned", {
  sim <- make_mid_rotation_sim()
  mv <- inia_merch_vol(sim, age = 16, pruned_height = 1e3)
  totals <- attr(mv, "totals")
  for (nm in c("vol_large_sawlog", "vol_small_sawlog", "vol_pulp")) {
    expect_equal(totals[[paste0(nm, "_unpruned")]], 0, tolerance = 1e-6)
    expect_equal(totals[[paste0(nm, "_pruned")]],
                 totals[[nm]], tolerance = 1e-6)
  }
})

test_that("pruned_height = 0 pushes all volume to unpruned", {
  sim <- make_mid_rotation_sim()
  mv <- inia_merch_vol(sim, age = 16, pruned_height = 0)
  totals <- attr(mv, "totals")
  for (nm in c("vol_large_sawlog", "vol_small_sawlog", "vol_pulp")) {
    expect_equal(totals[[paste0(nm, "_pruned")]], 0, tolerance = 1e-6)
    expect_equal(totals[[paste0(nm, "_unpruned")]],
                 totals[[nm]], tolerance = 1e-6)
  }
})

test_that("pruned_height rejects negative values", {
  sim <- make_mid_rotation_sim()
  expect_error(inia_merch_vol(sim, age = 16, pruned_height = -2),
               "non-negative")
})

test_that("compare_inia_grandis returns both sim results and a tidy frame", {
  cmp <- compare_inia_grandis(
    SI = 28, N0 = 900, G0 = 7,
    Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
    zone = 7, PASW = 130, Elev = 130,
    t0 = 2, t_end = 12, label = "site1"
  )
  expect_named(cmp, c("comparison", "inia", "grandis"))
  expect_true(is.data.frame(cmp$comparison))
  expect_true(all(c("INIA", "augmented") %in% cmp$comparison$model))
  expect_true("label" %in% names(cmp$comparison))
  expect_true(all(cmp$comparison$label == "site1"))
})

test_that("comparison frame is balanced: same number of rows per model", {
  cmp <- compare_inia_grandis(SI = 28, N0 = 900, G0 = 7,
                              Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
                              zone = 7, PASW = 130, Elev = 130,
                              t0 = 2, t_end = 14)
  n_per_model <- table(cmp$comparison$model)
  expect_equal(unname(n_per_model[1]), unname(n_per_model[2]))
})

test_that("comparison frame omits the label column when label is NULL", {
  cmp <- compare_inia_grandis(SI = 28, N0 = 900, G0 = 7,
                              Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
                              zone = 7, PASW = 130, Elev = 130,
                              t0 = 2, t_end = 10)
  expect_false("label" %in% names(cmp$comparison))
})

test_that("comparison stacks correctly across multiple sites", {
  sites <- data.frame(
    label = c("low", "high"),
    PASW  = c(80, 170),
    Elev  = c(200, 80)
  )
  combined <- do.call(rbind, lapply(seq_len(nrow(sites)), function(i) {
    compare_inia_grandis(
      SI = 28, N0 = 900, G0 = 7, Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
      zone = 7, PASW = sites$PASW[i], Elev = sites$Elev[i],
      t0 = 2, t_end = 12, label = sites$label[i]
    )$comparison
  }))
  # 4 series (2 sites x 2 models) on 11 ages
  expect_equal(nrow(combined), 4 * 11)
  expect_setequal(unique(combined$label), c("low", "high"))
})

test_that("thinning regimes are applied to both simulators", {
  cmp <- compare_inia_grandis(
    SI = 28, N0 = 900, G0 = 7, Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
    zone = 7, PASW = 130, Elev = 130,
    t0 = 2, t_end = 14,
    thins = list(list(age = 6, N_after = 500))
  )
  expect_equal(nrow(cmp$inia$thinnings), 1)
  expect_equal(nrow(cmp$grandis$thinnings), 1)
  expect_equal(cmp$inia$thinnings$N_post, 500)
  expect_equal(cmp$grandis$thinnings$N_post, 500)
})

test_that("plot_inia_grandis_compare returns a ggplot with the expected layers", {
  skip_if_not_installed("ggplot2")
  cmp <- compare_inia_grandis(
    SI = 28, N0 = 900, G0 = 7, Hd0 = 7, dmax0 = 13, SDd0 = 1.8,
    zone = 7, PASW = 130, Elev = 130, t0 = 2, t_end = 12, label = "demo"
  )
  p <- plot_inia_grandis_compare(cmp$comparison)
  expect_s3_class(p, "ggplot")
  # Two geom layers (line + point) plus the facet scales
  expect_gte(length(p$layers), 2)
})

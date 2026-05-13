# ------------------------------------------------------------------------------
# Diameter distribution: inverse Weibull recovered from stand moments
# Kuru et al. (1992); Methol (2001, 2003).
# ------------------------------------------------------------------------------


# Recover inverse Weibull parameters (a, b, c) from Dq, dmax, SDd via the
# method-of-moments approach in Methol 2001 / SAG 2003 eqn 17-18. Internal.
inia_weibull_params <- function(Dq, dmax, SDd) {
  Dmean <- sqrt(max(Dq^2 - SDd^2, 0.01))

  a <- dmax
  z_val <- SDd / max(a - Dmean, 0.01)
  z_val <- min(max(z_val, 0.01), 0.99)

  kz <- .inia_params$weibull_kz
  poly_sum <- sum(kz * z_val^(0:5))
  c_param <- 1 / (z_val * (1 + (1 - z_val)^2 * poly_sum))
  c_param <- max(c_param, 1.01)

  b_param <- (a - Dmean) / gamma(1 + 1 / c_param)
  b_param <- max(b_param, 0.01)

  list(a = a, b = b_param, c = c_param)
}

# Inverse Weibull CDF: F(DAP) = exp(-((a - DAP)/b)^c) for 0 <= DAP < a.
inia_weibull_cdf <- function(DAP, a, b, c) {
  if (DAP >= a) return(1)
  if (DAP <= 0) return(0)
  exp(-((a - DAP) / b)^c)
}


#' Diameter distribution from stand-level moments
#'
#' Recovers an inverse Weibull distribution from quadratic mean diameter,
#' maximum diameter, and diameter SD, then returns per-class frequencies
#' and basal area per hectare.
#'
#' @param N Trees per hectare.
#' @param Dq Quadratic mean diameter (cm).
#' @param dmax Maximum diameter (cm) - Weibull location parameter.
#' @param SDd Standard deviation of diameters (cm).
#' @param Hd Dominant height (m). Currently unused; reserved for future
#'   per-class height/volume estimation.
#' @param class_width Diameter class width (cm). Default 1.
#' @return Data frame with columns `class_mid` (cm), `freq` (trees/ha),
#'   `ba_tree` (m2 per tree), and `ba_ha` (m2/ha). Classes with negligible
#'   frequency (< 0.05 trees/ha) are dropped.
#' @examples
#' inia_diam_dist(N = 700, Dq = 28.0, dmax = 35.0, SDd = 4.5)
#' @export
inia_diam_dist <- function(N, Dq, dmax, SDd, Hd = NULL, class_width = 1) {
  params <- inia_weibull_params(Dq, dmax, SDd)
  a <- params$a; b <- params$b; c <- params$c

  classes <- seq(class_width, ceiling(dmax) + 5, by = class_width)

  result <- data.frame(class_mid = classes)

  result$freq <- vapply(classes, function(cl) {
    upper <- cl + class_width / 2
    lower <- cl - class_width / 2
    F_upper <- inia_weibull_cdf(upper, a, b, c)
    F_lower <- inia_weibull_cdf(lower, a, b, c)
    max((F_upper - F_lower) * N, 0)
  }, numeric(1))

  result$ba_tree <- pi * (result$class_mid / 200)^2
  result$ba_ha <- result$freq * result$ba_tree

  result[result$freq > 0.05, , drop = FALSE]
}

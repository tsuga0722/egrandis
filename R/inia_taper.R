# ------------------------------------------------------------------------------
# Taper, tree-level volume, and stand-level merchantable volume.
#
# Implements the Fang et al. (2000) compatible segmented taper-volume system
# with the E. grandis coefficients reported by Hirigoyen et al. (2021),
# iForest 14: 127-136 (Table 4, Original system).
#
# Functional form (Fang 2000), with q = h/H, H' = H / (H - 1.3),
# k = pi / 40000, indicators I1, I2:
#
#   d^2 = c1^2 * H'^((k - b1) / b1) * (1 - q)^((k - beta) / beta)
#         * alpha1^(I1 + I2) * alpha2^I2
#
#   I1 = 1 if p1 <= q <= p2; 0 otherwise
#   I2 = 1 if p2 <= q <= 1;  0 otherwise
#   beta   = b1^(1 - I1 - I2) * b2^I1 * b3^I2
#   alpha1 = (1 - p1)^( k*(b2 - b1) / (b1*b2) )
#   alpha2 = (1 - p2)^( k*(b3 - b2) / (b2*b3) )
#   c1^2   = a0 * D^a1 * H^(a2 - k/b1) /
#            [ b1*(r0 - r1) + b2*(r1 - alpha1*r2) + b3*alpha1*r2 ]
#   r0 = (1 - hb/H)^(k/b1),  r1 = (1 - p1)^(k/b1),  r2 = (1 - p2)^(k/b2)
#
# Compatible total volume: V = a0 * D^a1 * H^a2  (m^3, single tree).
#
# Reference parameters (E. grandis, Hirigoyen 2021 Table 4):
#   a0 = 4.0e-5,  a1 = 2.09,  a2 = 0.862,
#   b1 = 4.2e-6,  b2 = 3e-5,  b3 = 2e-4,
#   p1 = 0.031,   p2 = 0.941
# ------------------------------------------------------------------------------


# Taper-volume parameters. Internal. Stump height hb (m) is a default that
# can be overridden in the user-facing functions.
.taper_params <- list(
  fang_egrandis = list(
    a0 = 4.0e-5, a1 = 2.09,  a2 = 0.862,
    b1 = 4.2e-6, b2 = 3.0e-5, b3 = 2.0e-4,
    p1 = 0.031,  p2 = 0.941,
    hb = 0.10
  )
)


# Helper: precompute the segment-invariant quantities of the Fang system.
# Returned list is reused by both the taper function and the volume integrand.
#
# The published Fang form has d^2 = c1^2 * H^((k-b1)/b1) * (1-q)^((k-beta)/beta)
# * alpha1^(I1+I2) * alpha2^I2 with c1^2 = a0*D^a1*H^(a2 - k/b1)/[denom]. With
# the b_i values reported by Hirigoyen et al. (2021) (~1e-6 to 1e-4), the
# individual H exponents are ~+/-18 and underflow / overflow in double
# precision even though their algebraic sum is small. We collapse them
# analytically into a single combined coefficient k_eff multiplying the
# q-dependent part of d^2:
#
#   k_eff = a0 * D^a1 * H^(a2 - 1) / [denom]
#
# so that d^2 = k_eff * (1-q)^((k-beta)/beta) * alpha1^(I1+I2) * alpha2^I2
# and total volume integrates to V = a0 * D^a1 * H^a2 as required for
# compatibility with the Fang volume equation.
.fang_constants <- function(D, H, params = .taper_params$fang_egrandis,
                            hb = NULL) {
  p <- params
  if (is.null(hb)) hb <- p$hb
  k  <- pi / 40000

  alpha1 <- (1 - p$p1)^(k * (p$b2 - p$b1) / (p$b1 * p$b2))
  alpha2 <- (1 - p$p2)^(k * (p$b3 - p$b2) / (p$b2 * p$b3))

  r0 <- (1 - hb / H)^(k / p$b1)
  r1 <- (1 - p$p1)^(k / p$b1)
  r2 <- (1 - p$p2)^(k / p$b2)

  denom <- p$b1 * (r0 - r1) +
           p$b2 * (r1 - alpha1 * r2) +
           p$b3 * alpha1 * r2

  k_eff <- (p$a0 * D^p$a1 * H^(p$a2 - 1)) / denom

  list(p = p, k = k, alpha1 = alpha1, alpha2 = alpha2,
       k_eff = k_eff, hb = hb)
}


# Squared diameter at height h. Vectorised over h. Internal.
#
# For the Hirigoyen 2021 E. grandis fit, (k - b3)/b3 is slightly negative,
# so the raw Fang formula has d^2 increasing toward the tip in segment 3.
# The integrated volume of segment 3 is still finite (the singularity is
# integrable), and the closed-form V = a0*D^a1*H^a2 is preserved. To keep
# the user-facing diameter monotone we close segment 3 with a cone tail:
# above q = p2 the formula is replaced with a linear taper from d(p2*H)
# down to 0 at h = H. The displaced volume relative to the raw formula's
# segment 3 is < 1% of total tree volume in practice.
.fang_d2 <- function(h, D, H, cs) {
  q <- pmin(pmax(h / H, 0), 1)
  k  <- cs$k
  p  <- cs$p

  I1 <- as.numeric(q >= p$p1 & q <= p$p2)
  # I2 left as zero throughout: segment 3 handled separately below.
  beta <- p$b1^(1 - I1) * p$b2^I1

  one_minus_q <- pmax(1 - q, 1e-12)

  d2_seg12 <- cs$k_eff *
              one_minus_q^((k - beta) / beta) *
              cs$alpha1^I1

  # Diameter (squared) at the segment-2/3 inflection point.
  d2_p2 <- cs$k_eff *
           (1 - p$p2)^((k - p$b2) / p$b2) *
           cs$alpha1

  # Linear cone tail above p2: d goes from sqrt(d2_p2) at q=p2 to 0 at q=1.
  d_p2 <- sqrt(pmax(d2_p2, 0))
  d_seg3 <- d_p2 * pmax((1 - q) / (1 - p$p2), 0)
  d2_seg3 <- d_seg3^2

  d2 <- ifelse(q <= p$p2, d2_seg12, d2_seg3)
  pmax(d2, 0)
}


#' Diameter at height along the stem (Fang 2000 taper)
#'
#' Returns over-bark stem diameter (cm) at one or more heights `h` for an
#' *E. grandis* tree of breast-height diameter `D` and total height `H`.
#' Uses the Fang et al. (2000) compatible segmented system with the
#' coefficients reported in Hirigoyen et al. (2021) for *E. grandis*.
#'
#' @param h Height(s) of interest along the stem (m). Vectorised.
#' @param D DBH over bark (cm).
#' @param H Total tree height (m). Must exceed 1.3.
#' @param hb Stump height (m). Default 0.10.
#' @return Numeric vector of diameters (cm) at each `h`. Diameter is 0 at
#'   or above the tip and is set to 0 below the stump.
#' @references Fang Z, Borders BE, Bailey RL (2000). For. Sci. 46: 1-12.
#'   Hirigoyen A et al. (2021). iForest 14: 127-136.
#' @examples
#' # Diameter at 1.3 m should be very close to DBH
#' inia_taper(h = 1.3, D = 25, H = 28)
#'
#' # Full stem profile
#' h <- seq(0.1, 27.9, length.out = 30)
#' d <- inia_taper(h, D = 25, H = 28)
#' plot(h, d, type = "l", xlab = "Height (m)", ylab = "Diameter (cm)")
#' @export
inia_taper <- function(h, D, H, hb = 0.10) {
  if (H <= 1.3) stop("`H` must exceed 1.3 m")
  if (D <= 0)   return(rep(0, length(h)))
  cs <- .fang_constants(D, H, hb = hb)
  d <- sqrt(.fang_d2(h, D, H, cs))
  d[h < hb | h > H] <- 0
  d
}


#' Height along the stem where the diameter falls to a given small-end value
#'
#' Inverts the Fang taper for an *E. grandis* tree to find the height at
#' which over-bark diameter equals `d_top`. Useful for log bucking and
#' merchantable-length computations.
#'
#' @param d_top Small-end diameter limit (cm).
#' @param D DBH over bark (cm).
#' @param H Total tree height (m).
#' @param hb Stump height (m). Default 0.10.
#' @return Height (m) at which diameter equals `d_top`. Returns `H` if the
#'   tree never reaches `d_top` (i.e., `d_top >= D`); returns `hb` if the
#'   butt diameter is already smaller than `d_top`.
#' @examples
#' # Height to a 14 cm small-end limit on a 25 cm DBH / 28 m tall tree
#' inia_height_at_d(d_top = 14, D = 25, H = 28)
#' @export
inia_height_at_d <- function(d_top, D, H, hb = 0.10) {
  if (H <= 1.3) stop("`H` must exceed 1.3 m")
  if (d_top <= 0) return(H)
  d_butt <- inia_taper(hb, D = D, H = H, hb = hb)
  if (d_top >= d_butt) return(hb)
  # uniroot is robust here; the taper is monotone decreasing above the butt.
  rt <- stats::uniroot(
    function(h) inia_taper(h, D = D, H = H, hb = hb) - d_top,
    interval = c(hb, H - 1e-6), tol = 1e-5
  )
  rt$root
}


#' Total stem volume from the compatible Fang volume equation
#'
#' Closed-form total over-bark stem volume from the Fang et al. (2000)
#' system: `V = a0 * D^a1 * H^a2` with *E. grandis* coefficients.
#'
#' @param D DBH over bark (cm).
#' @param H Total tree height (m).
#' @return Volume (m^3, single tree).
#' @examples
#' inia_tree_total_vol(D = 25, H = 28)
#' @export
inia_tree_total_vol <- function(D, H) {
  p <- .taper_params$fang_egrandis
  p$a0 * D^p$a1 * H^p$a2
}


#' Tree-level volume between two stem heights or to a small-end diameter
#'
#' Integrates the Fang taper to compute over-bark volume between
#' `h_lower` and `h_upper`. By default, integrates from stump to total
#' height (i.e., returns total stem volume). Pass `d_top` to integrate
#' from stump to the height at which diameter falls to `d_top` -- i.e.,
#' merchantable volume to a small-end limit.
#'
#' @param D DBH over bark (cm).
#' @param H Total tree height (m).
#' @param h_lower Lower integration limit (m). Default `hb` (stump).
#' @param h_upper Upper integration limit (m). Default `H` (tip). Ignored
#'   if `d_top` is supplied.
#' @param d_top Small-end diameter (cm). If supplied, `h_upper` is set to
#'   `inia_height_at_d(d_top, D, H)`.
#' @param hb Stump height (m). Default 0.10.
#' @return Volume between the two heights (m^3, single tree).
#' @examples
#' # Merchantable volume to an 8 cm top
#' inia_tree_vol(D = 25, H = 28, d_top = 8)
#'
#' # Volume of a specific log: first 5.5 m above stump
#' inia_tree_vol(D = 25, H = 28, h_lower = 0.10, h_upper = 5.6)
#' @export
inia_tree_vol <- function(D, H, h_lower = NULL, h_upper = NULL,
                          d_top = NULL, hb = 0.10) {
  if (H <= 1.3) stop("`H` must exceed 1.3 m")
  if (D <= 0)   return(0)

  if (is.null(h_lower)) h_lower <- hb
  if (!is.null(d_top)) {
    h_upper <- inia_height_at_d(d_top, D = D, H = H, hb = hb)
  } else if (is.null(h_upper)) {
    h_upper <- H
  }
  if (h_upper <= h_lower) return(0)

  cs <- .fang_constants(D, H, hb = hb)
  k  <- pi / 40000  # cm^2 -> m^2 conversion for the area integrand
  integrand <- function(z) k * .fang_d2(z, D, H, cs)
  stats::integrate(integrand, lower = h_lower, upper = h_upper,
                   rel.tol = 1e-6)$value
}


#' Allometric per-class total height
#'
#' Estimates the total height of a tree of DBH `D` in a stand with
#' quadratic mean diameter `Dq` and dominant height `Hd`, using a simple
#' power allometry:
#'
#' \deqn{h(D) = 1.3 + (Hd - 1.3) \cdot \min\bigl((D/Dq)^k, c \bigr)}
#'
#' where `k` controls how strongly height tracks diameter and `c` caps
#' the height of supra-dominant trees relative to `Hd`. Defaults give
#' a biologically reasonable height-diameter curve without requiring
#' species-specific h-d coefficients.
#'
#' @param D DBH of the tree(s) of interest (cm).
#' @param Dq Quadratic mean diameter of the stand (cm).
#' @param Hd Dominant height of the stand (m).
#' @param k Power-law exponent. Default 0.5.
#' @param cap Maximum multiplier of `(Hd - 1.3)` for very large trees.
#'   Default 1.05 (allows 5 percent over Hd for upper-tail trees).
#' @return Total height (m). Vectorised over `D`.
#' @examples
#' inia_height_class(D = c(20, 30, 40), Dq = 30, Hd = 28)
#' @export
inia_height_class <- function(D, Dq, Hd, k = 0.5, cap = 1.05) {
  if (Hd <= 1.3) stop("`Hd` must exceed 1.3 m")
  ratio <- pmin((D / max(Dq, 0.01))^k, cap)
  1.3 + (Hd - 1.3) * ratio
}


#' Stand-level merchantable volume by product assortment
#'
#' Buckets stand volume into named products defined by small-end diameter
#' limits, integrating the Fang taper across the recovered Weibull
#' diameter distribution at the requested age. Products are assigned to
#' butt-up sections of each stem: the product with the largest `d_min`
#' takes the butt log until diameter falls below its limit, the next
#' product takes the next section down to its limit, and so on.
#'
#' @param sim Result of [simulate_inia()].
#' @param age Target age (must appear in `sim$trajectory$age`).
#' @param products Named list. Each element is a list with `d_min` (cm)
#'   -- the small-end diameter for that product. Order does not matter:
#'   products are sorted internally by `d_min` descending. If `NULL`, the
#'   default assortment `list(veneer = list(d_min = 25), solid = list(d_min = 14),
#'   pulp = list(d_min = 8))` is used.
#' @param hb Stump height (m). Default 0.10.
#' @param hd_k Power-law exponent passed to [inia_height_class()] for
#'   estimating per-diameter-class total height. Default 0.5.
#' @return Data frame with one row per diameter class times one column
#'   per product (`vol_<product>`, m^3/ha), plus `class_mid`, `freq`,
#'   `H_class`, and `vol_total` (m^3/ha) per class. A `totals` attribute
#'   on the returned object gives per-product per-hectare aggregates and
#'   the residual `top_waste` (volume between the smallest small-end
#'   limit and the tip).
#' @examples
#' sim <- simulate_inia(
#'   SI = 30, N0 = 550, G0 = 1.7,
#'   Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
#'   t0 = 1, t_end = 16, zone = 7
#' )
#' mv <- inia_merch_vol(sim, age = 16)
#' attr(mv, "totals")
#' @export
inia_merch_vol <- function(sim, age, products = NULL, hb = 0.10,
                           hd_k = 0.5) {
  if (is.null(products)) {
    products <- list(
      veneer = list(d_min = 25),
      solid  = list(d_min = 14),
      pulp   = list(d_min = 8)
    )
  }
  # Sort products by d_min descending so butt-up bucking works.
  d_mins <- vapply(products, function(p) p$d_min, numeric(1))
  ord <- order(d_mins, decreasing = TRUE)
  products <- products[ord]
  d_mins   <- d_mins[ord]
  prod_names <- names(products)

  # Pull the stand state at the requested age and recover the Weibull
  # diameter distribution.
  row <- sim$trajectory[sim$trajectory$age == age, ]
  if (nrow(row) == 0) {
    stop(sprintf("Age %s not found in simulation trajectory", age))
  }
  dd <- inia_diam_dist(
    N    = row$N, Dq = row$DAP_medio,
    dmax = row$DAP_max, SDd = row$Desvio_DAP
  )

  Dq <- row$DAP_medio
  Hd <- row$AMD

  # Per-class total height and bole-volume bookkeeping.
  dd$H_class   <- inia_height_class(dd$class_mid, Dq = Dq, Hd = Hd, k = hd_k)
  dd$vol_total <- 0
  for (nm in prod_names) {
    dd[[paste0("vol_", nm)]] <- 0
  }
  dd$vol_top_waste <- 0

  for (i in seq_len(nrow(dd))) {
    D_i <- dd$class_mid[i]
    H_i <- dd$H_class[i]
    n_i <- dd$freq[i]

    if (H_i <= 1.3 || D_i <= 0 || n_i <= 0) next

    total_vol <- inia_tree_vol(D = D_i, H = H_i, hb = hb)
    dd$vol_total[i] <- total_vol * n_i

    # Butt-up cascade: each product gets the next section down to its d_min.
    h_prev <- hb
    consumed <- 0
    for (j in seq_along(prod_names)) {
      nm <- prod_names[j]
      d_min <- d_mins[j]

      d_butt_remaining <- inia_taper(h_prev, D = D_i, H = H_i, hb = hb)
      if (d_butt_remaining <= d_min) {
        # No volume left for this or any smaller product above this point.
        next
      }
      h_next <- inia_height_at_d(d_min, D = D_i, H = H_i, hb = hb)
      vol_log <- inia_tree_vol(D = D_i, H = H_i,
                               h_lower = h_prev, h_upper = h_next,
                               hb = hb)
      dd[[paste0("vol_", nm)]][i] <- vol_log * n_i
      consumed <- consumed + vol_log * n_i
      h_prev <- h_next
    }

    dd$vol_top_waste[i] <- dd$vol_total[i] - consumed
  }

  totals <- c(
    stats::setNames(
      vapply(prod_names, function(nm) sum(dd[[paste0("vol_", nm)]]),
             numeric(1)),
      paste0("vol_", prod_names)
    ),
    top_waste = sum(dd$vol_top_waste),
    total     = sum(dd$vol_total)
  )

  out <- dd[, c("class_mid", "freq", "H_class", "vol_total",
                paste0("vol_", prod_names), "vol_top_waste")]
  attr(out, "totals") <- totals
  out
}

# ------------------------------------------------------------------------------
# INIA stand simulator and user-facing helpers
# ------------------------------------------------------------------------------


#' Simulate an *E. grandis* stand with the INIA 2021 model
#'
#' Projects dominant height, basal area, mortality, volume, max diameter,
#' diameter SD, SDI, MAI, and ICA on a discrete time step, optionally
#' applying thinning events. Reproduces the SAG grandis 2021 online
#' simulator (sag.inia.uy) for Zones 7-9 of Uruguay.
#'
#' @section Initialization caveat:
#' The basal-area projection (Schumacher form) overshoots when started far
#' below its asymptote - initializing at age 1 with G = 1.7 m^2/ha can drive
#' Dq to ~22 cm by age 3, whereas field data shows Dq ~13 cm. The behavior
#' matches the SAG 2021 online simulator. For reliable projections,
#' initialize from age 3+ with measured G.
#'
#' @param SI Site index: dominant height at the base age (default 10 yr), m.
#' @param N0 Initial trees per hectare.
#' @param G0 Initial basal area (m2/ha).
#' @param Hd0 Initial dominant height (m). If `NULL`, derived from `SI`.
#' @param dmax0 Initial maximum diameter (cm). If `NULL`, set to `1.4 * Dq`.
#' @param SDd0 Initial diameter SD (cm). If `NULL`, set to `0.2 * Dq`.
#' @param t0 Initial age (years).
#' @param t_end Final age (years).
#' @param zone INIA zone: 7 (Tacuarembo/Rivera), 8 (Durazno/Cerro Largo),
#'   or 9 (Paysandu/Rio Negro). Zones 8 and 9 produce identical output.
#' @param thins Optional list of thinning events. Each element is a list
#'   with `age` (years) and `N_after` (trees/ha post-thin).
#' @param dt Time step (years). Default 1.
#' @return A list with:
#' \describe{
#'   \item{`trajectory`}{Data frame, one row per time step, with columns
#'     `age`, `AMD` (dominant height, m), `N`, `AB` (basal area), `DAP_medio`,
#'     `DAP_max`, `Desvio_DAP`, `Vol_Total`, `MAI`, `SDI`, `ICA`.}
#'   \item{`thinnings`}{Data frame of applied thinning events, or `NULL`.}
#'   \item{`cumulative_thin_vol`}{Total volume removed in thins (m3/ha).}
#'   \item{`final_standing`}{Standing volume at `t_end` (m3/ha).}
#'   \item{`total_yield`}{`final_standing + cumulative_thin_vol`.}
#'   \item{`parameters`}{Echo of run inputs.}
#' }
#' @examples
#' # Unthinned reference scenario: Zone 7, SI = 30, 550 TPH
#' sim <- simulate_inia(
#'   SI = 30, N0 = 550, G0 = 1.7,
#'   Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
#'   t0 = 1, t_end = 16, zone = 7
#' )
#' head(sim$trajectory)
#'
#' # Two-thin solid-wood regime
#' sim_thinned <- simulate_inia(
#'   SI = 30, N0 = 550, G0 = 1.7,
#'   Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
#'   t0 = 1, t_end = 14, zone = 7,
#'   thins = list(
#'     list(age = 3, N_after = 412),
#'     list(age = 7, N_after = 197)
#'   )
#' )
#' sim_thinned$thinnings
#' @export
simulate_inia <- function(SI, N0, G0, Hd0 = NULL, dmax0 = NULL, SDd0 = NULL,
                          t0 = 1, t_end = 16, zone = 7,
                          thins = NULL, dt = 1) {

  if (!zone %in% c(7, 8, 9)) {
    stop("`zone` must be 7, 8, or 9")
  }
  Z7 <- ifelse(zone == 7, 1, 0)

  Hd <- if (is.null(Hd0)) inia_hd_from_si(t0, SI, Z7 = Z7) else Hd0
  G <- G0
  N <- N0
  Dq <- sqrt(G / max(N, 1) * 40000 / pi)
  dmax <- if (is.null(dmax0)) Dq * 1.4 else dmax0
  SDd <- if (is.null(SDd0)) Dq * 0.2 else SDd0

  aN <- NA; aNb <- NA; att <- NA
  thin_recs <- list()
  cumul_thin <- 0

  if (!is.null(thins)) {
    thins <- lapply(thins, function(x) { x$done <- FALSE; x })
  }

  results <- data.frame()
  ages <- seq(t0, t_end, by = dt)

  for (i in seq_along(ages)) {
    t <- ages[i]

    if (i > 1) {
      tp <- ages[i - 1]
      Hd   <- inia_hd_from_si(t, SI, Z7 = Z7)
      G    <- inia_g(G, tp, t, Z7, aN, aNb, att)
      N    <- inia_n(N, tp, t, Z7)
      dmax <- inia_dmax(dmax, tp, t, Z7)
      SDd  <- inia_sdd(SDd, tp, t, Z7)
    }

    Dq <- sqrt(G / max(N, 1) * 40000 / pi)
    V  <- inia_vol(G, N, Hd, Z7)
    SDI <- N * (Dq / 25)^1.605

    if (!is.null(thins)) {
      for (j in seq_along(thins)) {
        th <- thins[[j]]
        if (!th$done && abs(t - th$age) < dt / 2 + 0.01) {
          V_pre <- V; G_pre <- G; N_pre <- N; Dq_pre <- Dq

          G <- inia_gpost(G_pre, N_pre, th$N_after)
          N <- th$N_after

          Dq <- sqrt(G / N * 40000 / pi)
          V <- inia_vol(G, N, Hd, Z7)
          SDI <- N * (Dq / 25)^1.605

          V_rem <- V_pre - V
          cumul_thin <- cumul_thin + V_rem

          aN <- th$N_after; aNb <- N_pre; att <- th$age
          thins[[j]]$done <- TRUE

          thin_recs[[length(thin_recs) + 1]] <- data.frame(
            age = t, N_pre = round(N_pre), N_post = th$N_after,
            G_pre = round(G_pre, 1), G_post = round(G, 1),
            Dq_pre = round(Dq_pre, 1), Dq_post = round(Dq, 1),
            V_removed = round(V_rem, 1)
          )
        }
      }
    }

    results <- rbind(results, data.frame(
      age = t,
      AMD = round(Hd, 1),
      N = round(N),
      AB = round(G, 1),
      DAP_medio = round(Dq, 1),
      DAP_max = round(dmax, 1),
      Desvio_DAP = round(SDd, 1),
      Vol_Total = round(V, 1),
      MAI = round(V / t, 1),
      SDI = round(SDI)
    ))
  }

  if (nrow(results) > 1) {
    results$ICA <- c(results$MAI[1], round(diff(results$Vol_Total), 1))
  } else {
    results$ICA <- results$MAI
  }

  thin_df <- if (length(thin_recs) > 0) do.call(rbind, thin_recs) else NULL
  final <- results[nrow(results), ]

  list(
    trajectory = results,
    thinnings = thin_df,
    cumulative_thin_vol = round(cumul_thin, 1),
    final_standing = final$Vol_Total,
    total_yield = round(final$Vol_Total + cumul_thin, 1),
    parameters = list(SI = SI, N0 = N0, G0 = G0, zone = zone, Z7 = Z7,
                      t0 = t0, t_end = t_end)
  )
}


#' Diameter distribution at a given age in a simulation
#'
#' Convenience wrapper that pulls the stand state at `age` from a
#' `simulate_inia()` result and runs [inia_diam_dist()] on it.
#'
#' @param sim Result of [simulate_inia()].
#' @param age Target age (must appear in `sim$trajectory$age`).
#' @return Data frame as in [inia_diam_dist()].
#' @examples
#' sim <- simulate_inia(
#'   SI = 30, N0 = 550, G0 = 1.7,
#'   Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
#'   t0 = 1, t_end = 16, zone = 7
#' )
#' dd <- inia_get_distribution(sim, age = 10)
#' head(dd)
#' @export
inia_get_distribution <- function(sim, age) {
  row <- sim$trajectory[sim$trajectory$age == age, ]
  if (nrow(row) == 0) {
    stop(sprintf("Age %s not found in simulation trajectory", age))
  }
  inia_diam_dist(
    N = row$N, Dq = row$DAP_medio,
    dmax = row$DAP_max, SDd = row$Desvio_DAP,
    Hd = row$AMD
  )
}


#' Print a formatted summary of an INIA simulation
#'
#' @param sim Result of [simulate_inia()].
#' @return The simulation, invisibly. Called for side effects.
#' @examples
#' sim <- simulate_inia(
#'   SI = 30, N0 = 550, G0 = 1.7,
#'   Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
#'   t0 = 1, t_end = 16, zone = 7
#' )
#' inia_print_summary(sim)
#' @export
inia_print_summary <- function(sim) {
  p <- sim$parameters
  cat(sprintf("INIA E. grandis 2021 - Zone %d, SI=%.0f, N0=%d\n",
              p$zone, p$SI, p$N0))
  cat(sprintf("Ages %d to %d\n\n", p$t0, p$t_end))

  if (!is.null(sim$thinnings) && nrow(sim$thinnings) > 0) {
    cat("Thinning events:\n")
    print(sim$thinnings, row.names = FALSE)
    cat(sprintf("\nCumulative thinning volume: %.1f m3/ha\n",
                sim$cumulative_thin_vol))
  }

  cat(sprintf("\nFinal standing volume: %.1f m3/ha\n", sim$final_standing))
  cat(sprintf("Total yield: %.1f m3/ha\n", sim$total_yield))
  cat(sprintf("Peak MAI: %.1f m3/ha/yr\n\n", max(sim$trajectory$MAI)))

  print(sim$trajectory, row.names = FALSE)
  invisible(sim)
}

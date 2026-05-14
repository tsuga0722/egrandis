# ------------------------------------------------------------------------------
# Augmented simulator engine for E. grandis (`simulate_grandis()`)
#
# Orchestrates the literature-faithful submodels in grandis_model.R into
# a stand projection compatible with the INIA trajectory schema. See
# CLAUDE.md "Augmented Simulator" section for the design.
# ------------------------------------------------------------------------------


#' Simulate an *E. grandis* stand with the augmented simulator
#'
#' Projects dominant height, basal area, mortality (with optional
#' Reineke self-thinning), volume, max diameter, diameter SD, MAI, and
#' ICA on a discrete time step using literature-faithful submodels
#' calibrated on Uruguay PSPs (Rachid-Casnati et al. 2019 augmented
#' equations for Hd, G, dmax, SDd; Methol 2003 Eqn 6 for background
#' mortality; Reineke / Drew & Flewelling soft-logistic for
#' self-thinning; Fang taper x recovered Weibull for volume). All site
#' variation enters through the continuous site-variable inputs
#' (`PASW`, `Elev`, `slope`, `aspect`) rather than a zone dummy.
#'
#' @section Companion to `simulate_inia()`:
#' The trajectory schema is identical to [simulate_inia()] so the
#' downstream tools ([inia_merch_vol()], [inia_add_biomass()],
#' [inia_dmd_plot()], [inia_get_distribution()]) work unchanged. The
#' two simulators diverge in three places: the augmented equations
#' allow site-specific predictions where SAG 2021 / INIA has a zone
#' dummy; mortality has Reineke self-thinning superimposed on the
#' background; and volume is integrated from the Fang taper rather
#' than fit at the stand level.
#'
#' @section Provenance:
#' The returned `sim$provenance` list documents the per-submodel
#' source, the disclosed approximations, and an echo of the run
#' inputs. See `vignette("getting-started")` for guidance on choosing
#' between the two simulators.
#'
#' @param SI Site index: dominant height at the base age (default
#'   `base_age = 10` yr), m. RC2019 calibration envelope ~22-35.
#' @param N0 Initial trees per hectare.
#' @param G0 Initial basal area (m^2/ha).
#' @param PASW Potentially available soil water at the site (mm).
#'   RC2019 calibration envelope roughly 50-170 mm.
#' @param Elev Elevation of the site (m). RC2019 envelope ~50-250 m.
#' @param slope Slope of the site (percent). Default 0 (flat).
#' @param aspect Aspect (azimuth) of the slope, radians from north
#'   (0 = N, pi/2 = E, pi = S, 3*pi/2 = W). Default 0.
#' @param Hd0 Initial dominant height (m). If `NULL`, derived from
#'   `SI`, `base_age`, and the site variables.
#' @param dmax0 Initial maximum diameter (cm). If `NULL`, set to
#'   `1.4 * Dq` from `G0` and `N0`.
#' @param SDd0 Initial diameter SD (cm). If `NULL`, set to `0.2 * Dq`.
#' @param t0 Initial age (years).
#' @param t_end Final age (years). Default 20 (longer than the INIA
#'   default of 16 to cover veneer rotations).
#' @param thins Optional list of thinning events. Each element is a
#'   list with `age` (years) and `N_after` (trees/ha post-thin). Same
#'   shape as [simulate_inia()].
#' @param dt Time step (years). Default 1.
#' @param exo_rate Annual exogenous-mortality rate (fraction of trees
#'   lost per year to non-competition causes: wind, frost, drought,
#'   pests, lightning). Default `NULL` (uses `.grandis_params$mortality$exo_rate`
#'   = 0.005, i.e. 0.5%/year).
#' @param RD_ceiling Reineke threshold for the soft-logistic
#'   self-thinning blend. Default `NULL` (uses `.grandis_params$reineke$RD_50`
#'   = 0.60).
#' @param mort_k Slope of the soft-logistic mortality intensity in
#'   `RD - RD_ceiling`. Default `NULL` (uses 12; steep transition).
#' @param bark_factor Volume bark factor applied to the over-bark
#'   Fang integration. Default `NULL` (uses 0.82).
#' @param base_age Reference age for site index (years). Default 10.
#' @return A list with:
#' \describe{
#'   \item{`trajectory`}{Data frame with one row per time step and the
#'     same columns as [simulate_inia()]: `age`, `AMD`, `N`, `AB`,
#'     `DAP_medio`, `DAP_max`, `Desvio_DAP`, `Vol_Total`, `MAI`,
#'     `SDI`, `RD`, `ICA`.}
#'   \item{`thinnings`}{Data frame of applied thinning events, or `NULL`.}
#'   \item{`cumulative_thin_vol`}{Total volume removed in thins (m^3/ha).}
#'   \item{`final_standing`}{Standing volume at `t_end` (m^3/ha).}
#'   \item{`total_yield`}{`final_standing + cumulative_thin_vol`.}
#'   \item{`parameters`}{Echo of run inputs.}
#'   \item{`provenance`}{Named list documenting per-submodel sources,
#'     disclosed approximations, and run inputs. See Provenance section.}
#' }
#' @examples
#' # Default site (PASW=120, Elev=130, flat ground), SI=28, N0=900
#' sim <- simulate_grandis(
#'   SI = 28, N0 = 900, G0 = 7,
#'   PASW = 120, Elev = 130,
#'   t0 = 2, t_end = 18
#' )
#' head(sim$trajectory)
#' sim$provenance$caveats
#'
#' # Higher-quality site (more PASW, lower elevation, NE-facing) with
#' # a single commercial thin at age 8
#' sim_better <- simulate_grandis(
#'   SI = 30, N0 = 900, G0 = 7,
#'   PASW = 160, Elev = 100, slope = 6, aspect = pi/4,
#'   t0 = 2, t_end = 20,
#'   thins = list(list(age = 8, N_after = 450))
#' )
#' sim_better$total_yield
#' @export
simulate_grandis <- function(SI, N0, G0,
                             PASW, Elev,
                             slope = 0, aspect = 0,
                             Hd0 = NULL, dmax0 = NULL, SDd0 = NULL,
                             t0 = 1, t_end = 20,
                             thins = NULL,
                             dt = 1,
                             exo_rate = NULL,
                             RD_ceiling = NULL,
                             mort_k = NULL,
                             bark_factor = NULL,
                             base_age = 10) {

  # Resolve defaults pulled from the params table
  rp <- .grandis_params$reineke
  if (is.null(exo_rate))    exo_rate    <- .grandis_params$mortality$exo_rate
  if (is.null(RD_ceiling))  RD_ceiling  <- rp$RD_50
  if (is.null(mort_k))      mort_k      <- rp$mort_k
  if (is.null(bark_factor)) bark_factor <- .grandis_params$volume$bark_factor

  # Light input validation
  if (PASW < 0)   stop("`PASW` must be non-negative.")
  if (Elev < 0)   stop("`Elev` must be non-negative.")
  if (slope < 0)  stop("`slope` must be non-negative (in percent).")

  # Initial state
  Hd  <- if (is.null(Hd0)) grandis_hd_from_si(t0, SI, base_age,
                                              PASW = PASW,
                                              slope = slope, aspect = aspect) else Hd0
  G   <- G0
  N   <- N0
  Dq  <- sqrt(G / max(N, 1) * 40000 / pi)
  dmax <- if (is.null(dmax0)) Dq * 1.4 else dmax0
  SDd  <- if (is.null(SDd0)) Dq * 0.2 else SDd0

  # Thinning history slots for the G-modifier (latest event only)
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
      Hd   <- grandis_hd_from_si(t, SI, base_age,
                                  PASW = PASW, slope = slope, aspect = aspect)
      G    <- grandis_g(G, tp, t,
                        PASW = PASW, Elev = Elev,
                        Na = aN, Nb = aNb, tt = att)
      # Mortality consumes the projected G to compute Dq for the
      # Reineke check; project N AFTER G.
      N    <- grandis_n(N, tp, t, G = G,
                        exo_rate = exo_rate,
                        RD_50 = RD_ceiling, mort_k = mort_k)
      dmax <- grandis_dmax(dmax, tp, t,
                           PASW = PASW, slope = slope, aspect = aspect)
      SDd  <- grandis_sdd(SDd, tp, t, slope = slope, aspect = aspect)
    }

    Dq  <- sqrt(G / max(N, 1) * 40000 / pi)
    V   <- grandis_vol(G, N, Hd, dmax, SDd, bark_factor = bark_factor)
    SDI <- inia_sdi(N, Dq)

    if (!is.null(thins)) {
      for (j in seq_along(thins)) {
        th <- thins[[j]]
        if (!th$done && abs(t - th$age) < dt / 2 + 0.01) {
          V_pre <- V; G_pre <- G; N_pre <- N; Dq_pre <- Dq

          # Reuse INIA's Gpost (Methol 2003 Eqn 9, exact match)
          G  <- inia_gpost(G_pre, N_pre, th$N_after)
          N  <- th$N_after
          Dq <- sqrt(G / N * 40000 / pi)
          V  <- grandis_vol(G, N, Hd, dmax, SDd, bark_factor = bark_factor)
          SDI <- inia_sdi(N, Dq)

          V_rem <- V_pre - V
          cumul_thin <- cumul_thin + V_rem

          aN  <- th$N_after; aNb <- N_pre; att <- th$age
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
      SDI = round(SDI),
      RD = round(SDI / .inia_params$sdi$SDImax, 3)
    ))
  }

  results$ICA <- if (nrow(results) > 1) {
    c(results$MAI[1], round(diff(results$Vol_Total), 1))
  } else {
    results$MAI
  }

  thin_df <- if (length(thin_recs) > 0) do.call(rbind, thin_recs) else NULL
  final <- results[nrow(results), ]

  list(
    trajectory          = results,
    thinnings           = thin_df,
    cumulative_thin_vol = round(cumul_thin, 1),
    final_standing      = final$Vol_Total,
    total_yield         = round(final$Vol_Total + cumul_thin, 1),
    parameters          = list(
      SI = SI, N0 = N0, G0 = G0,
      PASW = PASW, Elev = Elev, slope = slope, aspect = aspect,
      t0 = t0, t_end = t_end, dt = dt, base_age = base_age,
      exo_rate = exo_rate,
      RD_ceiling = RD_ceiling, mort_k = mort_k, bark_factor = bark_factor
    ),
    provenance = .grandis_provenance(
      SI = SI, N0 = N0, G0 = G0,
      PASW = PASW, Elev = Elev, slope = slope, aspect = aspect,
      exo_rate = exo_rate,
      RD_ceiling = RD_ceiling, mort_k = mort_k, bark_factor = bark_factor,
      base_age = base_age
    )
  )
}


# Build the per-call provenance object. Internal.
.grandis_provenance <- function(...) {
  args <- list(...)
  list(
    model    = "simulate_grandis",
    version  = utils::packageVersion("egrandis"),
    submodels = list(
      Hd        = "RC2019 Eqn 11 (Tab S4 row 11) -- PASW + alpha_c + alpha_s; Uruguay PSPs",
      G         = "RC2019 Eqn 12 (Tab S5 row 12) -- PASW + Elev; Uruguay PSPs",
      dmax      = "RC2019 Eqn 13 (Tab S6 row 13) -- PASW + alpha_s; Uruguay PSPs",
      SDd       = "RC2019 Eqn 14 (Tab S7 row 14) -- alpha_s; Uruguay PSPs",
      N_bg      = "Constant-rate exogenous mortality (default 0.5%/year); user-overridable via `exo_rate`",
      N_self_thin = "Reineke (1933) + Drew & Flewelling (1979) soft-logistic; SDImax=1250 (Rachid-Casnati et al. 2024)",
      V         = "Fang taper (Hirigoyen 2021 iForest) integrated across recovered Weibull (Methol 2003)",
      Gpost     = "Methol 2003 Eqn 9 (empirical, reused from simulate_inia)",
      Weibull   = "Methol 2003 Eqns 10-11, 16-18 (reused from simulate_inia)"
    ),
    caveats = c(
      "All growth submodels calibrated on Uruguay PSPs (Rachid-Casnati et al. 2019, 305 plots, ages 2-11). Predictions outside that envelope are extrapolations.",
      "Mortality has two parts: a constant-rate exogenous term (default 0.5%/year, density-independent) and a Reineke / Drew-Flewelling self-thinning ceiling. Both rates are user-overridable but not species-fitted.",
      "Bark factor 0.82 is a published mean for E. grandis; varies with tree size and can be overridden.",
      "No local PSP calibration has been applied. simulate_grandis() is designed for future refit when PSP data is available."
    ),
    inputs = args
  )
}

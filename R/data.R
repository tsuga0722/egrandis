#' SAG grandis 2021 validation scenarios
#'
#' Output captured from the SAG grandis 2021 online simulator
#' (sag.inia.uy) in May 2026 for five scenarios spanning variation in
#' site index, initial density, INIA zone, and thinning regime. Used by
#' the `egrandis` test suite to validate that `simulate_inia()` reproduces
#' the reference simulator within documented tolerances.
#'
#' @format A named list with five elements. Each element is itself a list:
#' \describe{
#'   \item{`description`}{Human-readable description of the scenario.}
#'   \item{`inputs`}{Named list of input arguments suitable for passing to
#'     `simulate_inia()`.}
#'   \item{`trajectory`}{Data frame, one row per simulated age, with columns
#'     `age`, `AMD`, `N`, `AB`, `DAP_medio`, `DAP_max`, `Desvio_DAP`,
#'     `Vol_Total`, `MAI`, `ICA`, `SDI`.}
#'   \item{`thinnings`}{Data frame of thinning events, or `NULL` if
#'     unthinned.}
#' }
#'
#' The five scenarios are:
#' \itemize{
#'   \item `z7_si30_n550` - Zone 7, SI=30, 550 TPH, unthinned (reference).
#'   \item `z7_si25_n550` - Zone 7, SI=25, 550 TPH, unthinned.
#'   \item `z7_si30_n1111` - Zone 7, SI=30, 1111 TPH, unthinned (pulp).
#'   \item `z8_si30_n550` - Zone 8, SI=30, 550 TPH, unthinned.
#'   \item `z7_si30_n550_thinned` - Zone 7, SI=30, 550 TPH, two-thin
#'     solid-wood regime (age 3 -> 412, age 7 -> 197).
#' }
#'
#' @source SAG grandis 2021 online simulator, INIA Uruguay
#'   (\url{https://sag.inia.uy}), captured May 2026.
"sag_validation"

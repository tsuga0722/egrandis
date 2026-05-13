# Helpers used by the SAG validation tests.

# Run simulate_inia() against an entry of the sag_validation list,
# returning the matched trajectory rows from simulator and reference.
run_scenario <- function(name, skip_ages = integer(0)) {
  s <- sag_validation[[name]]
  args <- s$inputs
  args$DAP_medio0 <- NULL  # not a simulate_inia arg
  args$t0 <- 1
  args$t_end <- max(s$trajectory$age)
  sim <- do.call(simulate_inia, args)

  m <- sim$trajectory
  r <- s$trajectory
  ages <- setdiff(intersect(m$age, r$age), skip_ages)
  list(
    sim = sim,
    sim_traj = m[match(ages, m$age), , drop = FALSE],
    ref_traj = r[match(ages, r$age), , drop = FALSE],
    ages = ages
  )
}

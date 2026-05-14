# Generate the README figures into man/figures/ from package functions.
# Re-run after retuning parameters or any change that would update the
# README plots. Not run as part of R CMD check.

pkgload::load_all(".", quiet = TRUE)
out <- "man/figures"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

png_open <- function(name, w = 1200, h = 800) {
  grDevices::png(file.path(out, name), width = w, height = h, res = 150,
                 type = "cairo")
}


# --- Baseline scenario reused across figures -----------------------------------
baseline <- list(SI = 28, N0 = 900, G0 = 7.0,
                 Hd0 = 7.0, dmax0 = 13.0, SDd0 = 1.8,
                 t0 = 2, t_end = 16, zone = 7)
sim_baseline <- do.call(simulate_inia, baseline)


# --- 1. Taper profile of a typical mid-rotation tree ---------------------------
png_open("taper-profile.png", w = 1200, h = 800)
op <- par(mar = c(4.2, 4.2, 2.5, 1))
D <- 25; H <- 28
h_vals <- seq(0.1, H, length.out = 300)
d_vals <- inia_taper(h_vals, D = D, H = H)
plot(h_vals, d_vals, type = "l", lwd = 2, col = "darkgreen",
     xlab = "Height above ground (m)", ylab = "Stem diameter o.b. (cm)",
     main = sprintf("Fang taper, D = %d cm, H = %d m", D, H))
abline(h = c(8, 14, 25), lty = 3, col = "grey60")
text(rep(0.5, 3), c(8.7, 14.7, 25.7),
     c("8 cm (pulp)", "14 cm (small sawlog)", "25 cm (large sawlog)"),
     pos = 4, col = "grey30", cex = 0.9)
par(op); dev.off()


# --- 2. MAI curves by initial density (Zone 7, SI = 28) ------------------------
png_open("mai-by-density.png", w = 1200, h = 800)
op <- par(mar = c(4.2, 4.2, 2.5, 1))
densities <- c(800, 900, 1000, 1100)
g_from_dq <- function(N, Dq = 10) N * pi * (Dq / 200)^2
sims <- lapply(densities, function(n) {
  simulate_inia(SI = 28, N0 = n, G0 = g_from_dq(n),
                Hd0 = 7.0, dmax0 = 13.0, SDd0 = 1.8,
                t0 = 2, t_end = 16, zone = 7)
})
plot(NULL, xlim = c(2, 16), ylim = c(0, 50),
     xlab = "Age (years)", ylab = expression("MAI (m"^3*" ha"^-1*" yr"^-1*")"),
     main = "Mean annual increment by initial density (Zone 7, SI = 28)")
cols <- c("steelblue", "darkgreen", "darkorange", "firebrick")
for (i in seq_along(sims)) {
  lines(sims[[i]]$trajectory$age, sims[[i]]$trajectory$MAI,
        col = cols[i], lwd = 2)
}
legend("bottomright", legend = paste(densities, "TPH"),
       col = cols, lwd = 2, bty = "n")
par(op); dev.off()


# --- 3. Product assortment at age 16 -------------------------------------------
png_open("merch-assortment.png", w = 1200, h = 800)
op <- par(mar = c(5.5, 4.2, 2.5, 1))
mv <- inia_merch_vol(sim_baseline, age = 16)
t16 <- attr(mv, "totals")
bars <- c(`Large sawlog`  = unname(t16["vol_large_sawlog"]),
          `Small sawlog`  = unname(t16["vol_small_sawlog"]),
          `Pulp`          = unname(t16["vol_pulp"]),
          `Top waste`     = unname(t16["top_waste"]))
labels_under <- c(">=25 cm", "14-25 cm", "8-14 cm", "")
bp <- barplot(bars, col = c("#3B7BCE", "#E08821", "#5FA85A", "grey70"),
              ylab = expression("Volume (m"^3*" ha"^-1*")"),
              main = "Stand-level merch assortment at age 16",
              ylim = c(0, max(bars) * 1.18))
text(bp, bars + max(bars) * 0.04, sprintf("%.0f", bars), cex = 0.95)
mtext(labels_under, side = 1, at = bp, line = 2.2, cex = 0.85,
      col = "grey30")
par(op); dev.off()


# --- 4. Biomass / carbon / CO2 trajectory --------------------------------------
png_open("biomass-trajectory.png", w = 1200, h = 800)
sim_biom <- inia_add_biomass(sim_baseline)
op <- par(mar = c(4.2, 4.2, 2.5, 4.5))
with(sim_biom$trajectory, {
  plot(age, Biomasa, type = "l", lwd = 2, col = "darkgreen",
       xlab = "Age (years)", ylab = "Biomass (t/ha)",
       main = "Aboveground biomass and CO2 stock (Zone 7, SI = 28)",
       ylim = c(0, max(Biomasa) * 1.05))
  par(new = TRUE)
  plot(age, CO2eq, type = "l", lwd = 2, col = "steelblue", lty = 2,
       axes = FALSE, xlab = "", ylab = "")
  axis(4, col.axis = "steelblue", col = "steelblue")
  mtext("CO2 equivalent (t CO2/ha)", side = 4, line = 2.8,
        col = "steelblue")
  legend("topleft",
         legend = c("AGB (t/ha)", "CO2eq (t CO2/ha)"),
         col = c("darkgreen", "steelblue"), lwd = 2,
         lty = c(1, 2), bty = "n")
})
par(op); dev.off()


# --- 5. Density management diagrams (unthinned + thinned) --------------------
sim_thinned_dmd <- simulate_inia(
  SI = 28, N0 = 900, G0 = 7.0,
  Hd0 = 7.0, dmax0 = 13.0, SDd0 = 1.8,
  t0 = 2, t_end = 16, zone = 7,
  thins = list(
    list(age = 4, N_after = 600),
    list(age = 9, N_after = 300)
  )
)
ggplot2::ggsave(file.path(out, "dmd-unthinned.png"),
                inia_dmd_plot(sim_baseline),
                width = 8, height = 5.5, dpi = 150)
ggplot2::ggsave(file.path(out, "dmd-thinned.png"),
                inia_dmd_plot(sim_thinned_dmd),
                width = 8, height = 5.5, dpi = 150)


# --- 6. SAG validation overlay -------------------------------------------------
png_open("sag-validation.png", w = 1200, h = 800)
op <- par(mar = c(4.2, 4.2, 2.5, 1))
data(sag_validation, package = "egrandis")
ref <- sag_validation$z7_si30_n550
args <- ref$inputs; args$DAP_medio0 <- NULL
args$t0 <- 1; args$t_end <- 16
sim_v <- do.call(simulate_inia, args)
plot(ref$trajectory$age, ref$trajectory$Vol_Total,
     pch = 19, col = "black",
     xlab = "Age (years)", ylab = expression("Volume (m"^3*" ha"^-1*")"),
     main = "egrandis vs SAG 2021 reference scenario",
     ylim = range(c(ref$trajectory$Vol_Total, sim_v$trajectory$Vol_Total)))
lines(sim_v$trajectory$age, sim_v$trajectory$Vol_Total,
      col = "firebrick", lwd = 2)
legend("bottomright", c("SAG 2021 (reference)", "egrandis simulate_inia"),
       col = c("black", "firebrick"), pch = c(19, NA), lty = c(NA, 1),
       lwd = c(NA, 2), bty = "n")
par(op); dev.off()

cat("Generated:\n")
for (f in list.files(out, full.names = TRUE)) {
  cat("  ", f, " (",
      format(file.info(f)$size, big.mark = ","), " bytes)\n",
      sep = "")
}

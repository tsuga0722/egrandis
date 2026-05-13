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
     c("8 cm (pulp)", "14 cm (solid)", "25 cm (veneer)"),
     pos = 4, col = "grey30", cex = 0.9)
par(op); dev.off()


# --- 2. MAI curves by initial density (Zone 7, SI = 30) ------------------------
png_open("mai-by-density.png", w = 1200, h = 800)
op <- par(mar = c(4.2, 4.2, 2.5, 1))
densities <- c(400, 550, 800, 1111)
g_from_dq <- function(N, Dq = 6.3) N * pi * (Dq / 200)^2
sims <- lapply(densities, function(n) {
  simulate_inia(SI = 30, N0 = n, G0 = g_from_dq(n),
                Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
                t0 = 1, t_end = 16, zone = 7)
})
plot(NULL, xlim = c(1, 16), ylim = c(0, 55),
     xlab = "Age (years)", ylab = expression("MAI (m"^3*" ha"^-1*" yr"^-1*")"),
     main = "Mean annual increment by initial density (Zone 7, SI = 30)")
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
op <- par(mar = c(4.2, 4.2, 2.5, 1))
sim <- simulate_inia(SI = 30, N0 = 550, G0 = 1.7,
                     Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
                     t0 = 1, t_end = 16, zone = 7)
mv <- inia_merch_vol(sim, age = 16)
t16 <- attr(mv, "totals")
bars <- c(`Veneer >=25cm` = unname(t16["vol_veneer"]),
          `Solid 14-25cm` = unname(t16["vol_solid"]),
          `Pulp 8-14cm`   = unname(t16["vol_pulp"]),
          `Top waste`     = unname(t16["top_waste"]))
bp <- barplot(bars, col = c("#3B7BCE", "#E08821", "#5FA85A", "grey70"),
              ylab = expression("Volume (m"^3*" ha"^-1*")"),
              main = "Stand-level merch assortment at age 16 (550 TPH)",
              ylim = c(0, max(bars) * 1.18))
text(bp, bars + max(bars) * 0.04, sprintf("%.0f", bars), cex = 0.95)
par(op); dev.off()


# --- 4. Biomass / carbon / CO2 trajectory --------------------------------------
png_open("biomass-trajectory.png", w = 1200, h = 800)
sim <- inia_add_biomass(sim)
op <- par(mar = c(4.2, 4.2, 2.5, 4.5))
with(sim$trajectory, {
  plot(age, Biomasa, type = "l", lwd = 2, col = "darkgreen",
       xlab = "Age (years)", ylab = "Biomass (t/ha)",
       main = "Aboveground biomass and CO2 stock (Zone 7, SI = 30, 550 TPH)",
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


# --- 5. SAG validation overlay -------------------------------------------------
png_open("sag-validation.png", w = 1200, h = 800)
op <- par(mar = c(4.2, 4.2, 2.5, 1))
data(sag_validation, package = "egrandis")
ref <- sag_validation$z7_si30_n550
args <- ref$inputs; args$DAP_medio0 <- NULL
args$t0 <- 1; args$t_end <- 16
sim <- do.call(simulate_inia, args)
plot(ref$trajectory$age, ref$trajectory$Vol_Total,
     pch = 19, col = "black",
     xlab = "Age (years)", ylab = expression("Volume (m"^3*" ha"^-1*")"),
     main = "egrandis vs SAG 2021 (Z7, SI=30, 550 TPH unthinned)",
     ylim = range(c(ref$trajectory$Vol_Total, sim$trajectory$Vol_Total)))
lines(sim$trajectory$age, sim$trajectory$Vol_Total,
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

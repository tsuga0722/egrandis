# egrandis

**Growth and Yield Models for *Eucalyptus grandis* Clonal Plantations**

An R package for simulating growth, mortality, thinning response, and yield in *E. grandis* plantations across a range of management regimes. Designed for operational forestry planning in subtropical South America (Paraguay, Uruguay, Argentina, southern Brazil).

## Overview

`egrandis` provides two complementary modeling systems:

1. **INIA model** (`simulate_inia()`) — A faithful reproduction of the SAG grandis 2021 system developed by INIA Uruguay (Rachid-Casnati, Hirigoyen, Varela). Validated against the online simulator at sag.inia.uy. Best suited as a reference benchmark for managed thinning regimes within its calibration domain (500–1100 TPH, thinned, SI 22–35, Zone 7 northern Uruguay).

2. **Paraguay model** (`simulate_py()`) — *In development.* A state-space whole-stand model with density-dependent mortality, drought response, and a causal Hd → RD → N → G → V architecture. Designed for the full design space including low-density veneer regimes (400 TPH), unthinned carbon-max plantations (1111 TPH), and P-deficient Alfisol sites with drought risk.

## Installation

```r
# Install from GitHub
devtools::install_github("tsuga0722/egrandis")
```

## Quick Start

```r
library(egrandis)

# Simulate an unthinned stand: Zone 7, SI=30, 550 TPH
result <- simulate_inia(
  SI = 30, N0 = 550, G0 = 1.7,
  Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
  t0 = 1, t_end = 16, zone = 7
)

# View the trajectory
result$trajectory

# Print a formatted summary
inia_print_summary(result)

# Get the diameter distribution at age 10
dd <- inia_get_distribution(result, 10)
plot(dd$class_mid, dd$freq, type = "h",
     xlab = "DBH (cm)", ylab = "Trees/ha",
     main = "Diameter distribution at age 10")
```

### Thinning Regimes

```r
# Two-thin solid wood regime
result <- simulate_inia(
  SI = 30, N0 = 550, G0 = 1.7,
  Hd0 = 5.2, dmax0 = 8.0, SDd0 = 1.3,
  t0 = 1, t_end = 14, zone = 7,
  thins = list(
    list(age = 3, N_after = 412),
    list(age = 7, N_after = 197)
  )
)

inia_print_summary(result)
```

### Comparing Scenarios

```r
# Density comparison
densities <- c(400, 550, 800, 1111)
results <- lapply(densities, function(n) {
  g0 <- n * pi * (6.3 / 200)^2  # BA from Dq=6.3cm
  simulate_inia(SI = 30, N0 = n, G0 = g0, t0 = 1, t_end = 16, zone = 7)
})

# Plot MAI curves
plot(NULL, xlim = c(1, 16), ylim = c(0, 50),
     xlab = "Age (years)", ylab = "MAI (m3/ha/yr)")
cols <- c("blue", "darkgreen", "orange", "red")
for (i in seq_along(results)) {
  lines(results[[i]]$trajectory$age, results[[i]]$trajectory$MAI,
        col = cols[i], lwd = 2)
}
legend("topright", paste(densities, "TPH"), col = cols, lwd = 2)
```

## Model Details

### INIA Model Submodels

| Submodel | Equation | Source | Validation |
|----------|----------|--------|------------|
| Dominant height (Hd) | Johnson-Schumacher, polymorphic | Rachid-Casnati et al. 2019, eqn 5 | Exact match |
| Basal area (G) | Schumacher(2) + thinning modifier | Rachid-Casnati et al. 2019, eqn 6 | Exact match |
| Mortality (N) | Clutter-Jones | Fitted to SAG 2021 output | ±12 trees |
| Volume (V) | Stand-level exponential | Fitted to SAG 2021 output | ±2.5% |
| Post-thin BA | Empirical | Methol 2003, eqn 9 | ~2% |
| Max diameter (dmax) | Schumacher projection | Fitted to SAG 2021 output | ±1.2 cm |
| Diameter SD (SDd) | von Bertalanffy-Richards | Rachid-Casnati et al. 2019, eqn 8 | Exact match |
| Diameter distribution | Inverse Weibull | Methol 2001/2003 | ±1-3 trees/class |

### Zones

- **Zone 7** (Tacuarembó, Rivera — northern Uruguay): Primary parameterization, most relevant for Paraguay. Higher growth rates, SI range 22–35.
- **Zones 8/9** (central/western Uruguay): Lower productivity, identical output in SAG 2021. Not recommended for Paraguay applications.

### Known Limitations

The INIA model has specific limitations that users should understand:

- **Mortality is not density-dependent in its parameters.** The Clutter-Jones functional form creates implicit density dependence through the state variable (higher N₁ → proportionally more mortality), but there is no explicit stocking feedback. A 150 TPH stand loses trees at the same proportional rate as a 1100 TPH stand. This is appropriate for the calibration domain but unreliable at very low or very high densities.

- **Basal area overshoots when initialized at age 1.** The Schumacher projection form converges aggressively toward the asymptote (~56 m²/ha for Zone 7) from a low starting BA. Initializing at age 1 with G=1.7 produces G=19.3 at age 3 (Dq=22.4 cm), which is biologically unrealistic — field data shows Dq≈13 cm at age 3. For reliable projections, initialize from age 3+ with measured G, or use the age-1 initialization with the understanding that early BA is overstated.

- **No competition-driven self-thinning.** Unthinned stands at high density will not exhibit density-dependent mortality or self-thinning dynamics. All mortality is a background process.

- **Volume is under-bark, stand-level.** The volume equation approximates the SAG 2021's diameter-distribution-integrated volume with a stand-level shortcut. Systematic bias is ≈2.5% for ages 3+.

## Package Data

The package includes a built-in dataset `sag_validation` containing output from the SAG grandis 2021 online simulator for five scenarios, used for validation and examples:

```r
data(sag_validation)
str(sag_validation)
```

## References

- Rachid-Casnati C, Mason E, Woollons R (2019). Using soil-based and physiographic variables to improve stand growth equations in Uruguayan forest plantations. *iForest* 12: 237–245. doi:10.3832/ifor2926-012

- Methol R (2003). SAG grandis: Sistema de Apoyo a la Gestión de Plantaciones de *Eucalyptus grandis*. INIA Serie Técnica 131.

- Methol R (2001). Comparisons of approaches to modelling tree taper, stand structure and stand dynamics in forest plantations. PhD thesis, University of Canterbury.

- Rachid-Casnati C, Resquin F, Hirigoyen A et al. (2024). Impact of thinning on the yield and quality of *Eucalyptus grandis* wood at harvest time in Uruguay. *Forests* 15: 810.

- Resquin F, Navarro-Cerrillo RM, Rachid-Casnati C et al. (2018). Allometry, growth and survival of three *Eucalyptus* species in high-density plantations in Uruguay. *Forests* 9: 745.

## License

MIT

## Contributing

Contributions welcome. Please open an issue first to discuss proposed changes. If you have PSP data from *E. grandis* plantations in the region and are interested in model calibration or validation, we would be glad to collaborate.

## Run this script from the package root to regenerate README figures.
## The results are already hardcoded in README.md/README.Rmd — this only
## regenerates the PNGs in man/figures/.

library(sensemakr)
source("R/bootmakr.R")
data(darfur)

f <- peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
       pastvoted + hhsize_darfur + female + village

cat("Generating figures (this takes a few minutes)...\n\n")

# kd sweep, simple bootstrap
out2 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 benchmark_covariates = "female",
                 kd = seq(0.5, 3, by = 0.5),
                 reps = 5000, seed = 42, progress = TRUE)

# Cluster bootstrap
out3 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 benchmark_covariates = "female",
                 kd = seq(0.5, 3, by = 0.5),
                 reps = 5000, seed = 42,
                 cluster = "village", progress = TRUE)

# Grouped benchmark
out4 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 gbenchmark_covariates = c("female", "pastvoted"),
                 kd = seq(0.5, 3, by = 0.5),
                 reps = 5000, seed = 42,
                 cluster = "village", progress = TRUE)

# Convergence
out5 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 benchmark_covariates = "female", kd = 1,
                 reps = 5000, seed = 42,
                 cluster = "village",
                 converge = list(minreps = 500, stepsize = 500, threshold = 3000),
                 progress = TRUE)

# Save figures
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

png("man/figures/README-kd-sweep-1.png", width = 900, height = 550, res = 150, pointsize = 8)
plot(out2, type = "kd_sweep")
dev.off()

png("man/figures/README-cluster-1.png", width = 900, height = 550, res = 150, pointsize = 8)
plot(out3, type = "kd_sweep")
dev.off()

png("man/figures/README-gbenchmark-1.png", width = 900, height = 550, res = 150, pointsize = 8)
plot(out4, type = "kd_sweep")
dev.off()

png("man/figures/README-convergence-1.png", width = 900, height = 900, res = 150, pointsize = 8)
plot(out5, type = "convergence")
dev.off()

cat("\nAll figures saved to man/figures/\n")

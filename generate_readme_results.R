## Run this script to generate all README output and figures.
## Copy-paste the console output and send it back so I can hardcode it.

library(sensemakr)
source("R/bootmakr.R")   # adjust path if needed
data(darfur)

f <- peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
       pastvoted + hhsize_darfur + female + village


# ==============================================================================
# 1. Quick start: single kd, simple bootstrap
# ==============================================================================
cat("==== EXAMPLE 1: QUICK START ====\n")
out1 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 benchmark_covariates = "female", kd = 1,
                 reps = 5000, seed = 42, progress = TRUE)
print(out1)


# ==============================================================================
# 2. kd sweep, simple bootstrap (heteroskedasticity)
# ==============================================================================
cat("\n\n==== EXAMPLE 2: KD SWEEP ====\n")
out2 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 benchmark_covariates = "female",
                 kd = seq(0.5, 3, by = 0.5),
                 reps = 5000, seed = 42, progress = TRUE)
print(out2)


# ==============================================================================
# 3. Cluster bootstrap on village
# ==============================================================================
cat("\n\n==== EXAMPLE 3: CLUSTER BOOTSTRAP ====\n")
out3 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 benchmark_covariates = "female",
                 kd = seq(0.5, 3, by = 0.5),
                 reps = 5000, seed = 42,
                 cluster = "village", progress = TRUE)
print(out3)


# ==============================================================================
# 4. Grouped benchmark (female + pastvoted), cluster bootstrap
# ==============================================================================
cat("\n\n==== EXAMPLE 4: GROUPED BENCHMARK ====\n")
out4 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 gbenchmark_covariates = c("female", "pastvoted"),
                 kd = seq(0.5, 3, by = 0.5),
                 reps = 5000, seed = 42,
                 cluster = "village", progress = TRUE)
print(out4)


# ==============================================================================
# 5. Convergence diagnostics, cluster bootstrap
# ==============================================================================
cat("\n\n==== EXAMPLE 5: CONVERGENCE ====\n")
out5 <- bootmakr(f, data = darfur, treat = "directlyharmed",
                 benchmark_covariates = "female", kd = 1,
                 reps = 5000, seed = 42,
                 cluster = "village",
                 converge = list(minreps = 500, stepsize = 500, threshold = 3000),
                 progress = TRUE)
print(out5)


# ==============================================================================
# 6. Bootstrap draws summary
# ==============================================================================
cat("\n\n==== EXAMPLE 6: DRAWS ====\n")
draws <- out3$boot_samples[, 1]
draws <- draws[is.finite(draws)]
cat("quantile(draws, c(0.025, 0.5, 0.975)):\n")
print(quantile(draws, c(0.025, 0.5, 0.975)))


# ==============================================================================
# Figures
# ==============================================================================
cat("\n\nGenerating figures...\n")

png("man/figures/README-kd-sweep-1.png", width = 700, height = 450, res = 150)
plot(out2, type = "kd_sweep")
dev.off()

png("man/figures/README-cluster-1.png", width = 700, height = 450, res = 150)
plot(out3, type = "kd_sweep")
dev.off()

png("man/figures/README-gbenchmark-1.png", width = 700, height = 450, res = 150)
plot(out4, type = "kd_sweep")
dev.off()

png("man/figures/README-convergence-1.png", width = 700, height = 700, res = 150)
plot(out5, type = "convergence")
dev.off()

cat("All figures saved to man/figures/\n")
cat("\n==== ALL DONE ====\n")

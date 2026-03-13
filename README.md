
# bootmakr

<!-- badges: start -->

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Bootstrap inference for sensitivity analysis under omitted variable
bias.**

`bootmakr` wraps the
[sensemakr](https://github.com/carloscinelli/sensemakr) package in a
bootstrap loop, producing bootstrap standard errors, percentile
confidence intervals, and *p*-values for the bias-adjusted treatment
effect. It is the companion R package to the Stata command of the same
name.

## Why bootstrap the sensitivity bounds?

`sensemakr` computes analytical adjusted estimates and confidence
intervals under a hypothetical confounder with a given strength. These
analytical intervals rely on asymptotic OLS standard errors, which may
not perform well with clustered data, small samples, or complex survey
designs. `bootmakr` replaces the analytical inference with a
nonparametric bootstrap — including cluster and stratified bootstrap —
so the CIs and *p*-values are robust to these complications.

## Installation

``` r
# Install from GitHub (once published):
# devtools::install_github("username/bootmakr")
```

## Quick start

``` r
library(bootmakr)

data(darfur, package = "sensemakr")

out <- bootmakr(
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar + pastvoted,
  data    = darfur,
  treat   = "directlyharmed",
  benchmark_covariates = "age",
  kd      = 1,
  reps    = 1000,
  seed    = 42,
  progress = FALSE
)
out
#> 
#> Call:
#> bootmakr(formula = peacefactor ~ directlyharmed + age + farmer_dar + 
#>     herder_dar + pastvoted, data = darfur, treat = "directlyharmed", 
#>     benchmark_covariates = "age", kd = 1, reps = 1000, seed = 42, 
#>     progress = FALSE)
#> 
#> Bootstrap sensitivity analysis (1,000 reps, n = 1,276)
#> Benchmark: age | kd = 1, ky = 1
#> 
#> Adjusted estimates (percentile 95% CI):
#>                Estimate Std. Err     2.5%    97.5% Pr(>|0|)   
#> directlyharmed 0.060236 0.019465 0.022831 0.097697    0.004 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
```

## Sweeping across benchmark strengths

Supply a vector of `kd` values to see how the adjusted effect changes as
the hypothetical confounder grows stronger. The `plot()` method produces
a coefficient plot with bootstrap CIs:

``` r
out_sweep <- bootmakr(
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar + pastvoted,
  data    = darfur,
  treat   = "directlyharmed",
  benchmark_covariates = "age",
  kd      = seq(0.5, 3, by = 0.5),
  reps    = 1000,
  seed    = 42,
  progress = FALSE
)
out_sweep
#> 
#> Call:
#> bootmakr(formula = peacefactor ~ directlyharmed + age + farmer_dar + 
#>     herder_dar + pastvoted, data = darfur, treat = "directlyharmed", 
#>     benchmark_covariates = "age", kd = seq(0.5, 3, by = 0.5), 
#>     reps = 1000, seed = 42, progress = FALSE)
#> 
#> Bootstrap sensitivity analysis (1,000 reps, n = 1,276)
#> Benchmark: age | kd = 0.5 1 1.5 2 2.5 3, ky = 0.5 1 1.5 2 2.5 3
#> 
#> Adjusted estimates (percentile 95% CI):
#>                         Estimate Std. Err     2.5%    97.5% Pr(>|0|)   
#> directlyharmed (kd=0.5) 0.060720 0.019445 0.023569 0.098267    0.004 **
#> directlyharmed (kd=1)   0.060236 0.019465 0.022831 0.097697    0.004 **
#> directlyharmed (kd=1.5) 0.059752 0.019504 0.022081 0.097408    0.004 **
#> directlyharmed (kd=2)   0.059268 0.019562 0.021542 0.097064    0.006 **
#> directlyharmed (kd=2.5) 0.058784 0.019639 0.021089 0.096775    0.004 **
#> directlyharmed (kd=3)   0.058299 0.019735 0.020823 0.096479    0.004 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
plot(out_sweep, type = "kd_sweep")
```

![](man/figures/README-kd-sweep-1.png)<!-- -->

## Cluster bootstrap

Pass a cluster identifier to `cluster` for a cluster-robust bootstrap —
the resampling is done at the cluster level:

``` r
# Simulated clustered data
set.seed(999)
G  <- rep(1:50, each = 20)
u  <- rnorm(50)[G]
x  <- 0.5 * u + rnorm(1000)
cv <- 0.3 * u + rnorm(1000)
y  <- 0.5 * x + 0.4 * cv + 0.8 * u + rnorm(1000)
sim <- data.frame(y = y, x = x, cv = cv, cid = G)

out_cl <- bootmakr(
  y ~ x + cv, data = sim, treat = "x",
  benchmark_covariates = "cv", kd = c(1, 2, 3),
  reps = 1000, seed = 77, cluster = "cid",
  progress = FALSE
)
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
#> Warning in ovb_partial_r2_bound.numeric(r2dxj.x = r2dxj.x[i], r2yxj.dx =
#> r2yxj.dx[i], : Implied bound on r2yz.dx greater than 1, try lower kd and/or ky.
#> Setting r2yz.dx to 1.
out_cl
#> 
#> Call:
#> bootmakr(formula = y ~ x + cv, data = sim, treat = "x", benchmark_covariates = "cv", 
#>     kd = c(1, 2, 3), reps = 1000, seed = 77, cluster = "cid", 
#>     progress = FALSE)
#> 
#> Bootstrap sensitivity analysis (1,000 reps, n = 1,000, 50 clusters)
#> Benchmark: cv | kd = 1 2 3, ky = 1 2 3
#> 
#> Adjusted estimates (percentile 95% CI):
#>          Estimate Std. Err     2.5%    97.5% Pr(>|0|)    
#> x (kd=1) 0.755700 0.052849 0.647354 0.853459        0 ***
#> x (kd=2) 0.728366 0.051354 0.616263 0.821577        0 ***
#> x (kd=3) 0.700961 0.057128 0.578025 0.804919        0 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
```

## Grouped benchmarks

When the benchmark for the hypothetical confounder should reflect the
*joint* explanatory power of several covariates, use
`gbenchmark_covariates`. Internally this computes the group partial R²
(via `sensemakr::group_partial_r2`) and applies the proper nonlinear
kd-scaling:

``` r
out_g <- bootmakr(
  lg_ceopay ~ owner_ceo + ceo_tenure + PA_nic3_med + ceo_edu_dummy +
    lg_sales + firm_age + promoters_pct + institutions_pct + year + nic,
  data    = ceo_data,
  treat   = "owner_ceo",
  gbenchmark_covariates = c("ceo_tenure", "ceo_edu_dummy"),
  kd      = seq(0.1, 0.5, by = 0.1),
  reps    = 10000,
  cluster = "firm_id",
  seed    = 912323
)
plot(out_g, type = "kd_sweep")
```

## Convergence diagnostics

Large-sample bootstrap inference depends on using enough replications.
Pass `converge = TRUE` (or a list with fine-grained control) to assess
whether SEs and *p*-values have stabilised:

``` r
out_conv <- bootmakr(
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar + pastvoted,
  data    = darfur,
  treat   = "directlyharmed",
  benchmark_covariates = "age",
  kd      = 1,
  reps    = 2000,
  seed    = 42,
  converge = list(minreps = 200, stepsize = 200, threshold = 1500),
  progress = FALSE
)
out_conv
#> 
#> Call:
#> bootmakr(formula = peacefactor ~ directlyharmed + age + farmer_dar + 
#>     herder_dar + pastvoted, data = darfur, treat = "directlyharmed", 
#>     benchmark_covariates = "age", kd = 1, reps = 2000, seed = 42, 
#>     progress = FALSE, converge = list(minreps = 200, stepsize = 200, 
#>         threshold = 1500))
#> 
#> Bootstrap sensitivity analysis (2,000 reps, n = 1,276)
#> Benchmark: age | kd = 1, ky = 1
#> 
#> Adjusted estimates (percentile 95% CI):
#>                Estimate Std. Err     2.5%    97.5% Pr(>|0|)   
#> directlyharmed 0.060236 0.019539 0.022839 0.098752    0.003 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
#> 
#> Convergence diagnostics (reps 200 to 2000 by 200, threshold = 1500):
#>              Mean  Range  CV % Range (>=1500) CV% (>=1500)
#> Std. error 0.0194 0.0011  2.13         0.0002         0.57
#> P-value    0.0022 0.0040 71.91         0.0008        14.24
#> (CV % = coefficient of variation: sd / mean * 100)
plot(out_conv, type = "convergence")
```

![](man/figures/README-convergence-1.png)<!-- -->

## Accessing the raw bootstrap draws

All bootstrap replicates are stored in the returned object, so there is
no need to re-run the analysis to inspect the distribution:

``` r
draws <- out$boot_samples[, 1]          # vector of B adjusted estimates
draws <- draws[is.finite(draws)]

# Percentiles, moments, etc.
quantile(draws, c(0.025, 0.5, 0.975))
#>       2.5%        50%      97.5% 
#> 0.02283058 0.05960335 0.09769660

# Or export for further analysis
# write.csv(data.frame(adjusted_estimate = draws), "boot_draws.csv")
```

## Key arguments

| Argument | Description |
|----|----|
| `formula`, `data`, `treat` | Standard OLS specification and treatment name |
| `benchmark_covariates` | Individual benchmark covariate(s) |
| `gbenchmark_covariates` | Grouped benchmark covariates (joint partial R²) |
| `kd`, `ky` | Benchmark strength multipliers (`ky` defaults to `kd`) |
| `reps`, `seed` | Number of bootstrap replications and random seed |
| `cluster`, `strata` | Cluster and/or strata identifiers |
| `alpha` | Significance level (default 0.05) |
| `converge` | `TRUE`, `FALSE`, or `list(minreps, stepsize, threshold)` |
| `progress` | Show a progress bar (default `TRUE`) |

## Methods

| Method | Description |
|----|----|
| `print(x)` | Coefficient table with bootstrap SEs, CIs, and *p*-values |
| `plot(x, type = "kd_sweep")` | Coefficient plot across kd values |
| `plot(x, type = "histogram")` | Bootstrap distribution histogram |
| `plot(x, type = "convergence")` | Three-panel convergence diagnostic plot |
| `plot(x)` | Auto-selects the most informative plot |

## References

Cinelli, C. and Hazlett, C. (2020). Making Sense of Sensitivity:
Extending Omitted Variable Bias. *Journal of the Royal Statistical
Society, Series B (Statistical Methodology)*, 82(1), 39–67.

## License

MIT

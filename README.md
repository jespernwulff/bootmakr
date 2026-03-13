
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

We use the Darfur data from `sensemakr`, with the full model
specification including village fixed effects and `female` as the
benchmark covariate:

``` r
library(bootmakr)

data(darfur, package = "sensemakr")

out <- bootmakr(
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female + village,
  data    = darfur,
  treat   = "directlyharmed",
  benchmark_covariates = "female",
  kd      = 1,
  reps    = 500,
  seed    = 42,
  progress = FALSE
)
out
#> 
#> Call:
#> bootmakr(formula = peacefactor ~ directlyharmed + age + farmer_dar + 
#>     herder_dar + pastvoted + hhsize_darfur + female + village, 
#>     data = darfur, treat = "directlyharmed", benchmark_covariates = "female", 
#>     kd = 1, reps = 500, seed = 42, progress = FALSE)
#> 
#> Bootstrap sensitivity analysis (500 reps, n = 1,276)
#> Benchmark: female | kd = 1, ky = 1
#> 
#> Adjusted estimates (percentile 95% CI):
#>                Estimate Std. Err     2.5%    97.5% Pr(>|0|)  
#> directlyharmed 0.075220 0.026589 0.014209 0.119091    0.012 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
```

## Sweeping across benchmark strengths

Supply a vector of `kd` values to see how the adjusted effect changes as
the hypothetical confounder grows stronger. The `plot()` method produces
a coefficient plot with bootstrap CIs. Here the simple (non-clustered)
bootstrap accounts for heteroskedasticity induced by the village fixed
effects:

``` r
out_sweep <- bootmakr(
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female + village,
  data    = darfur,
  treat   = "directlyharmed",
  benchmark_covariates = "female",
  kd      = 1:3,
  reps    = 500,
  seed    = 42,
  progress = FALSE
)
out_sweep
#> 
#> Call:
#> bootmakr(formula = peacefactor ~ directlyharmed + age + farmer_dar + 
#>     herder_dar + pastvoted + hhsize_darfur + female + village, 
#>     data = darfur, treat = "directlyharmed", benchmark_covariates = "female", 
#>     kd = 1:3, reps = 500, seed = 42, progress = FALSE)
#> 
#> Bootstrap sensitivity analysis (500 reps, n = 1,276)
#> Benchmark: female | kd = 1 2 3, ky = 1 2 3
#> 
#> Adjusted estimates (percentile 95% CI):
#>                       Estimate Std. Err      2.5%    97.5% Pr(>|0|)  
#> directlyharmed (kd=1) 0.075220 0.026589  0.014209 0.119091    0.012 *
#> directlyharmed (kd=2) 0.052915 0.031626 -0.025391 0.098804     0.18  
#> directlyharmed (kd=3) 0.030396 0.038848 -0.062536 0.083173    0.656  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
plot(out_sweep, type = "kd_sweep")
```

![](man/figures/README-kd-sweep-1.png)<!-- -->

## Cluster bootstrap

If observations are correlated within villages, pass
`cluster = "village"` for a cluster-robust bootstrap — the resampling is
done at the village level:

``` r
out_cl <- bootmakr(
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female + village,
  data    = darfur,
  treat   = "directlyharmed",
  benchmark_covariates = "female",
  kd      = 1:3,
  reps    = 500,
  seed    = 42,
  cluster = "village",
  progress = FALSE
)
out_cl
#> 
#> Call:
#> bootmakr(formula = peacefactor ~ directlyharmed + age + farmer_dar + 
#>     herder_dar + pastvoted + hhsize_darfur + female + village, 
#>     data = darfur, treat = "directlyharmed", benchmark_covariates = "female", 
#>     kd = 1:3, reps = 500, seed = 42, cluster = "village", progress = FALSE)
#> 
#> Bootstrap sensitivity analysis (500 reps, n = 1,276, 486 clusters)
#> Benchmark: female | kd = 1 2 3, ky = 1 2 3
#> 
#> Adjusted estimates (percentile 95% CI):
#>                       Estimate Std. Err      2.5%    97.5% Pr(>|0|)   
#> directlyharmed (kd=1) 0.075220 0.025932  0.022021 0.125235    0.004 **
#> directlyharmed (kd=2) 0.052915 0.033005 -0.015680 0.114647     0.12   
#> directlyharmed (kd=3) 0.030396 0.043038 -0.062151 0.110465    0.452   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
plot(out_cl, type = "kd_sweep")
```

![](man/figures/README-cluster-1.png)<!-- -->

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
  reps    = 500,
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
  peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female + village,
  data    = darfur,
  treat   = "directlyharmed",
  benchmark_covariates = "female",
  kd      = 1,
  reps    = 1000,
  seed    = 42,
  cluster = "village",
  converge = list(minreps = 100, stepsize = 100, threshold = 750),
  progress = FALSE
)
out_conv
#> 
#> Call:
#> bootmakr(formula = peacefactor ~ directlyharmed + age + farmer_dar + 
#>     herder_dar + pastvoted + hhsize_darfur + female + village, 
#>     data = darfur, treat = "directlyharmed", benchmark_covariates = "female", 
#>     kd = 1, reps = 1000, seed = 42, cluster = "village", progress = FALSE, 
#>     converge = list(minreps = 100, stepsize = 100, threshold = 750))
#> 
#> Bootstrap sensitivity analysis (1,000 reps, n = 1,276, 486 clusters)
#> Benchmark: female | kd = 1, ky = 1
#> 
#> Adjusted estimates (percentile 95% CI):
#>                Estimate Std. Err     2.5%    97.5% Pr(>|0|)   
#> directlyharmed 0.075220 0.025740 0.021085 0.125726    0.006 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)
#> 
#> Convergence diagnostics (reps 100 to 1000 by 100, threshold = 750):
#>              Mean  Range  CV % Range (>=750) CV% (>=750)
#> Std. error 0.0265 0.0020  2.81        0.0003        0.52
#> P-value    0.0025 0.0067 99.83        0.0042       44.27
#> (CV % = coefficient of variation: sd / mean * 100)
plot(out_conv, type = "convergence")
```

![](man/figures/README-convergence-1.png)<!-- -->

## Accessing the raw bootstrap draws

All bootstrap replicates are stored in the returned object, so there is
no need to re-run the analysis to inspect the distribution:

``` r
draws <- out_cl$boot_samples[, 1]      # vector of B adjusted estimates
draws <- draws[is.finite(draws)]

# Percentiles, moments, etc.
quantile(draws, c(0.025, 0.5, 0.975))
#>       2.5%        50%      97.5% 
#> 0.02202080 0.07589556 0.12523548

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

Cinelli, C., J. Ferwerda, and C. Hazlett (2024). sensemakr: Sensitivity
analysis tools for OLS in R and Stata. *Observational Studies*, 10(2),
93-127.

Lonati, S. and J. N. Wulff (2026). Why you should not use the ITCV with
robust standard errors (and what to do instead). *SSRN Working Paper*.

## License

MIT

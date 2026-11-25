# bootmakr (development version)

This version lives on the `experimental` branch. The released package on
`master` is 0.4.0 and is not changed while the paper is under review.

* New: a **bootstrap robustness value**. `print()` now ends with the
  robustness values of the treatment coefficient, and they are stored in
  `$robustness`: `sensemakr`'s two values (the strength of an omitted
  variable, equally associated with treatment and outcome, that brings the
  adjusted estimate to zero, `rv_q`, and the strength at which it is no
  longer significant at `alpha` under the conventional OLS standard error,
  `rv_qa`) and their bootstrap counterpart, `rv_qa_boot`: the smallest
  strength at which the bootstrap test in the table of adjusted estimates
  stops rejecting. Running `bootmakr()` again with the same seed and
  `r2dz.x` and `r2yz.dx` marginally above that value gives a p-value of at
  least `alpha`, and marginally below it a p-value below `alpha`. It needs
  no benchmark and no extra model fits, works with every
  resampling scheme including two-way clustering, and is available for
  `q = 1` and `reduce = TRUE`. It has the same target as the table, the
  sampling variability of the adjusted estimate at a given strength of the
  omitted variable, and is not a statement about the robust standard error
  of the regression that would include that variable. As in the table, the
  omitted variable pulls the coefficient of each replicate towards zero
  whatever its sign there, which raises the bootstrap value when many
  replicates have the opposite sign (a weak effect).
* New element `$boot_base`: the unadjusted treatment coefficient and its
  scale (standard error times the square root of the residual degrees of
  freedom) in every replicate. The adjusted estimate at any strength of the
  omitted variable is a closed form in these two numbers.
* Adjusted estimates, bootstrap draws and all other results are unchanged
  for a given seed.
* With `weights`, a benchmark given through `benchmark_covariates` is now
  measured with the weights on both sides. Until now only the outcome model
  used them: the benchmark's partial R² with the treatment came from an
  unweighted regression, as in `sensemakr` 0.1.5, while `gbenchmark_covariates`
  already used the weights, so the same covariate could give two different
  answers. Adjusted estimates, bootstrap results and the benchmark-strength
  block of a weighted analysis with `benchmark_covariates` therefore change.
  Results without weights, and weighted results with `gbenchmark_covariates`,
  do not. `$sensemakr_orig` is still `sensemakr`'s own object, so its
  bounds can differ from those in `$results` when weights are used.
* The documentation of `weights` now says that weights are treated as fixed:
  each unit is resampled together with its weight, and weights are not
  re-estimated within replications.

# bootmakr 0.4.0

* The benchmark-strength block now reports the **impact** of the assumed
  omitted variable: the product of its partial correlations with the outcome
  and with the treatment, both given the covariates only, which is the scale
  on which the impact threshold of a confounding variable (ITCV) is stated.
  The outcome correlation is recovered from `sensemakr`'s bound (which also
  conditions on the treatment) with the recursion formula for partial
  correlations, signed so that the omitted variable biases the estimate away
  from zero (towards zero with `reduce = FALSE`). The column replaces the two
  absolute partial correlations in the printed block; both remain in
  `$benchmark_strength$implied`, which also gains `r_yz.x` and `impact`, and
  `$results` gains `impact`. `$benchmark_strength$r_yd.x` holds the partial
  correlation of the outcome with the treatment used in the recursion. Like
  the rest of the block, the impact is computed from the data alone.
* New example data set `firms`: a simulated panel of 250 firms observed for
  20 years, a variant of the simulated example in Lonati and Wulff (2026),
  with an observed control `c` and an unobserved variable `q` of equal
  strength, so that "an omitted variable as strong as `c`" is a variable like
  `q` and every answer can be checked. The examples, the vignette and the
  README use it instead of `sensemakr`'s Darfur data; the Stata command ships
  the same data as `firms.dta`.
* Documentation: the paper is cited as Lonati and Wulff (2026), *Academy of
  Management Proceedings*, <https://doi.org/10.5465/AMPROC.2026.247bp>, and
  the bootstrap procedure is credited to Cinelli, Ferwerda and Hazlett (2024,
  Appendix C). New website article, *Why bootstrap?*, with the paper's
  argument against the ITCV and analytic `sensemakr` under robust standard
  errors, and a worked example of what happens when `kd` is impossibly large.

# bootmakr 0.3.0

* `print()` now reports the **benchmark strength** below the table of
  adjusted estimates: the partial R² of the benchmark with the treatment and
  with the outcome, the strength of the omitted variable that `sensemakr`
  derives from it at each `kd`, and the corresponding partial correlations.
  The block is stored in `$benchmark_strength`, and the implied `r2dz.x` and
  `r2yz.dx` are added to `$results`. These quantities are computed from the
  data alone, so they do not depend on the variance estimator and bootstrap
  results for a given seed are unchanged.
* `weights` now works. Any formula created outside `bootmakr()` used to stop
  with "object 'weight_vec' not found".
* Replication counts of 10,000 or more are printed as integers ("10,000 reps"
  rather than "1e+04 reps").
* The returned object and the `print()` and `plot()` methods are documented.
* New vignette, "Get started with bootmakr in R", and a website with worked
  examples in R and Stata: <https://jespernwulff.github.io/bootmakr/>.

# bootmakr 0.2.0

* Two-way cluster bootstrap with the Cameron, Gelbach and Miller (2011)
  subtractive variance estimator: pass two cluster variables,
  `cluster = c("firm", "year")`.

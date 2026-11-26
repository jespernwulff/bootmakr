# bootmakr (development version)

This version lives on the `experimental` branch. The released package on
`master` is 0.4.0 and is not changed while the paper is under review.

* New: **standard errors side by side**. After the robustness values
  `print()` shows, for each `kd`, the bootstrap standard error of the table
  next to three others: the bootstrap standard error with the benchmark's
  strength held at its value in the original fit, computed from the same
  replicates; `sensemakr`'s adjusted standard error; and the closed-form
  robust standard error of Wainstein and Hazlett (2026, Option 3), which is
  HC1 without clusters, CR1 with one cluster variable and, as an extension
  made here, two-way CR1. The last two refer to the regression that would
  include the omitted variable, the bootstrap to the sampling variability of
  the adjusted estimate, so the columns need not agree. Stored in
  `$se_compare`.
* New: the **ratio of the robust to the conventional standard error** of the
  treatment coefficient, and its change when each benchmark is taken out of
  the regression and restored (Acedański, 2026, Section 6.2). Printed last
  and stored in `$se_ratio`.
* New argument `fix_benchmark`. With `TRUE` the benchmark's partial R² are
  taken from the original fit in every replicate instead of being
  re-estimated, so the omitted variable has the same strength throughout,
  as in the bootstrap of Wainstein and Hazlett (2026). The table and the
  breakdown kd then leave out the uncertainty about the benchmark.
* New plot, `plot(x, type = "contour")`: contours of the limit of the
  bootstrap interval nearer zero over the two partial R² of the omitted
  variable, from the stored replicates, with a line where the bootstrap test
  stops rejecting. The strength is held fixed across replicates, so at a
  benchmark point the interval differs from the table unless
  `fix_benchmark = TRUE`.
* New argument `twoway_crit`. With `"t"` the two-way interval, p-value,
  breakdown kd and bootstrap robustness value use a t distribution with
  min(G, H) - 1 degrees of freedom instead of the normal. The default is
  unchanged.
* `weights` may be a function of the data. It is then called on every
  resample, so the weights are re-estimated within replications (Option 1 of
  Wainstein and Hazlett, 2026). Benchmarks are measured with those weights
  and not with their semi-weights.
* Fixed: with `bounds_row` other than 1, or with a manual R² supplied next
  to a benchmark, the columns `r2dz.x`, `r2yz.dx` and `impact` of `$results`
  described the first benchmark's bound and not the row the adjusted
  estimate refers to. Adjusted estimates and their inference were not
  affected. The released 0.4.0 has the same mislabelling.
* New: the **breakdown kd**. After the benchmark-strength block `print()`
  reports the smallest `kd` (with `ky = kd`) at which the estimate is no
  longer significant at `alpha`, for the benchmark the table refers to:
  once judged by `sensemakr`'s interval with the conventional standard
  error, and once by the bootstrap test of the table. Running `bootmakr()`
  again with the same seed and a `kd` marginally below the bootstrap value
  gives a p-value below `alpha`, marginally above it a p-value of at least
  `alpha`. It is exact for the replicates drawn and needs no extra model
  fits or a `kd` sweep: the benchmark's two partial R² are now stored for
  every replicate (`$boot_bench`), and `sensemakr`'s bound and adjusted
  estimate are closed forms in them. The largest `kd` the benchmark allows
  is printed too. Replicates beyond their own limit are left out, as in the
  table, and a note says how many; near the limit the bias grows without
  bound, so a breakdown kd close to it reflects the limit more than the
  estimate. Available for `reduce = TRUE` and for all resampling schemes;
  stored in `$breakdown`.
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
* The breakdown kd and the bootstrap robustness value are set by the tail
  of the bootstrap distribution, so they need 1,000 or more replications
  and still differ by a few percent between seeds. In both, the adjustment
  follows the sign of each replicate's coefficient, as in the table, which
  can lift the bootstrap value above the conventional one when the effect
  is weak.
* New element `$boot_base`: the unadjusted treatment coefficient and its
  scale (standard error times the square root of the residual degrees of
  freedom) in every replicate. The adjusted estimate at any strength of the
  omitted variable is a closed form in these two numbers.
* Adjusted estimates, bootstrap draws and all other results are unchanged
  for a given seed, apart from the three columns named in the fix above.
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

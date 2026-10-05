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

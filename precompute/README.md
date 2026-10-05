# Precomputed vignette, articles and home page

The vignette, the website articles and the home page show the output of
bootstrap runs with up to 40,000 replications, in R and in Stata. They are
therefore not rendered from scratch when the package or the site is built.
Each one has a source file ending in `.Rmd.orig`. Knitting the source once
produces a plain `.Rmd` (or, for the home page, `.md`) that already contains
the code, its output and its figures. `R CMD build` and `pkgdown` only format
that file.

**Edit the `.Rmd.orig` files, never the generated ones.**

| Source | Generated | Stata runs it quotes |
|---|---|---|
| `vignettes/bootmakr.Rmd.orig` | `vignettes/bootmakr.Rmd` | none |
| `vignettes/articles/stata.Rmd.orig` | `vignettes/articles/stata.Rmd` | `first-steps` |
| `vignettes/articles/ceo-pay.Rmd.orig` | `vignettes/articles/ceo-pay.Rmd` | `first-steps`, `ceo-pay-1`, `ceo-pay-2`, `ceo-pay-3` |
| `vignettes/articles/benchmark-strength.Rmd.orig` | `vignettes/articles/benchmark-strength.Rmd` | none |
| `vignettes/articles/clustering.Rmd.orig` | `vignettes/articles/clustering.Rmd` | `first-steps` |
| `vignettes/articles/stata-help.Rmd.orig` | `vignettes/articles/stata-help.Rmd` | `stata-help` |
| `pkgdown/index.Rmd.orig` | `pkgdown/index.md` | `first-steps`, `ceo-pay-1` |

`vignettes/articles/r-and-stata.Rmd` has no computations and is written by
hand.

## Rebuilding

1. **Stata.** Start each run file from `precompute/stata/`:

   ```
   cd precompute/stata
   stata -e do first-steps.do
   stata -e do ceo-pay-1.do
   stata -e do ceo-pay-2.do
   stata -e do ceo-pay-3.do
   stata -e do stata-help.do
   ```

   The runs use the installed `bootmakr`, or the one in the folder named by
   the environment variable `BOOTMAKR_STATA_DIR` (a checkout of
   <https://github.com/jespernwulff/bootmakr-stata>); `stata-help.do` needs
   that variable. The `ceo-pay` runs read
   `CEO pay project - Data for sharing FINAL.dta`, the data of Chen, Chittoor
   and Vissa (2021), from the folder named by `BOOTMAKR_CHEN_DIR`;
   `first-steps.do` downloads the same file from the authors' OSF
   repository, as the R sources do. Everything the articles quote is written to
   `precompute/stata/output/`: the log (`<run>.txt`), key results, timings
   and graphs. The three `ceo-pay` runs take 10 to 45 minutes each and can
   run side by side.

2. **R.** With the development version of the package installed, from the
   package root:

   ```
   Rscript precompute/precompile.R            # everything
   Rscript precompute/precompile.R ceo-pay    # one source
   ```

   Chunks with long runs are cached in `precompute/cache/` (not under
   version control), so only changed chunks are run again. Needs `knitr`,
   `haven` and `sandwich`, and an internet connection for the example data.
   Figures are embedded in the generated files; the image files that knitr
   leaves in `vignettes/figures/` and `vignettes/articles/figures/` are not
   needed afterwards, except that the home page takes its figure from the
   worked example.

3. **Site.** `pkgdown::build_site()` writes the site to `docs/`, which is
   what GitHub Pages serves.

## How Stata output gets into the articles

The run files mark every block that an article shows:

```stata
* >>> kd1
bootmakr lg_ceopay owner_ceo $controls, ///
    treat(owner_ceo) benchmark(ceo_tenure) cluster(co_code) ///
    reps(10000) seed(912323)
* <<<
```

`stata_show("ceo-pay-1", "kd1")` in an article (see `helpers.R`) prints the
commands of that block, read from the do-file, and below them what Stata
printed, read from the log. No output is typed or pasted by hand, and the
numbers quoted in the running text are computed from the same results.

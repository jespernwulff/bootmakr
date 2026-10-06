# Generates the example data set `firms`: data/firms.rda for the R package and
# firms.dta for the Stata command (written next to this repository, into the
# clone of jespernwulff/bootmakr-stata, if that folder exists).
#
# The data are a variant of the simulated example in Lonati and Wulff (2026,
# Table 1). As there, G = 250 firms are observed for T = 20 years, a
# firm-level variable q affects both x and y, and the within-firm components
# of x and of the disturbance are first-order autoregressive with persistence
# .90, so that standard errors must be clustered by firm. Two things differ:
# an observed firm-level control c enters x and y exactly as q does, so that
# the benchmark approach of sensemakr applies (q is "an omitted variable as
# strong as c"), and the coefficients of c and q are 0.5 instead of 1 and 0.7,
# which keeps the implied strength of the omitted variable well inside the
# range sensemakr can handle. The random draws of q, v and u are those of seed
# 116 of the paper's script (c is drawn after them).
#
#   q_g ~ N(0, 1), c_g ~ N(0, 1)                        firm-level
#   v_gt, u_gt: AR(1) within firm, persistence .9, stationary variance 1
#   x_gt = 0.5 q_g + 0.5 c_g + v_gt
#   y_gt = 0.25 x_gt + 0.5 q_g + 0.5 c_g + u_gt
#
# Run from the package root:  Rscript data-raw/firms.R

G  <- 250
T_ <- 20

ar1_panel <- function(rho) {
  as.vector(vapply(seq_len(G), function(g) {
    z <- numeric(T_)
    z[1] <- rnorm(1)
    for (t in 2:T_) z[t] <- rho * z[t - 1] + rnorm(1, 0, sqrt(1 - rho^2))
    z
  }, numeric(T_)))
}

set.seed(116)
firm <- rep(seq_len(G), each = T_)
q    <- rep(rnorm(G), each = T_)
v    <- ar1_panel(0.9)
u    <- ar1_panel(0.9)
c_   <- rep(rnorm(G), each = T_)
x    <- 0.5 * q + 0.5 * c_ + v
y    <- 0.25 * x + 0.5 * q + 0.5 * c_ + u

firms <- data.frame(
  firm = as.integer(firm),
  year = rep(seq_len(T_), times = G),
  y    = y,
  x    = x,
  c    = c_,
  q    = q
)

save(firms, file = "data/firms.rda", version = 2, compress = "xz")
cat("data/firms.rda:", nrow(firms), "rows\n")

stata_dir <- file.path("..", "bootmakr-stata")
if (dir.exists(stata_dir) && requireNamespace("haven", quietly = TRUE)) {
  d <- firms
  attr(d$firm, "label") <- "Firm identifier (cluster)"
  attr(d$year, "label") <- "Year within firm, 1 to 20"
  attr(d$y, "label")    <- "Outcome"
  attr(d$x, "label")    <- "Predictor of interest; true effect on y is 0.25"
  attr(d$c, "label")    <- "Observed firm-level control"
  attr(d$q, "label")    <- "Unobserved firm-level variable (omitted in the analysis)"
  attr(d, "label") <- "bootmakr example data: variant of Lonati & Wulff (2026), Table 1"
  haven::write_dta(d, file.path(stata_dir, "firms.dta"), version = 14)
  cat("wrote", file.path(stata_dir, "firms.dta"), "\n")
}

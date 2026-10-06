# The impact of the omitted variable (its product of partial correlations with
# the outcome and the treatment, given the covariates only) is recovered from
# sensemakr's bound with the recursion formula for partial correlations. For
# an observed variable the result must equal the directly computed product.

data("firms", package = "bootmakr")

partial_cor <- function(a, b, z, data) {
  ra <- residuals(lm(reformulate(z, a), data = data))
  rb <- residuals(lm(reformulate(z, b), data = data))
  cor(ra, rb)
}

test_that("manual strength equal to q's own gives q's own impact", {
  r_xq.c <- partial_cor("x", "q", "c", firms)
  r_yq.c <- partial_cor("y", "q", "c", firms)
  fitq   <- lm(y ~ x + c + q, data = firms)
  r2yz.dx <- unname(sensemakr::partial_r2(fitq, covariates = "q"))

  out <- bootmakr(y ~ x + c, data = firms, treat = "x",
                  r2dz.x = r_xq.c^2, r2yz.dx = r2yz.dx,
                  reps = 20, seed = 1, progress = FALSE)
  bs <- out$benchmark_strength
  expect_identical(bs$type, "manual")

  fit <- lm(y ~ x + c, data = firms)
  expect_equal(bs$r_yd.x,
               unname(sign(coef(fit)[["x"]]) * sqrt(sensemakr::partial_r2(fit, covariates = "x"))),
               tolerance = 1e-10)
  expect_equal(bs$r_yd.x, partial_cor("y", "x", "c", firms), tolerance = 1e-10)

  # the recursion recovers q's partial correlation with y given c alone,
  # and the impact is the product with its partial correlation with x
  expect_equal(bs$implied$r_yz.x, r_yq.c, tolerance = 1e-8)
  expect_equal(bs$implied$impact, r_xq.c * r_yq.c, tolerance = 1e-8)
  expect_equal(out$results$impact, bs$implied$impact)
})

test_that("benchmark path: impact is consistent with the stored columns and with reduce", {
  out <- bootmakr(y ~ x + c, data = firms, treat = "x",
                  benchmark_covariates = "c", kd = c(0.5, 1, 2),
                  reps = 20, seed = 1, progress = FALSE)
  imp <- out$benchmark_strength$implied
  r_yd <- out$benchmark_strength$r_yd.x
  expect_true(r_yd > 0)
  r_yz.x <- imp$r_yz.dx * sqrt((1 - r_yd^2) * (1 - imp$r2dz.x)) + r_yd * imp$r_dz.x
  expect_equal(imp$r_yz.x, r_yz.x, tolerance = 1e-12)
  expect_equal(imp$impact, r_yz.x * imp$r_dz.x, tolerance = 1e-12)
  expect_equal(out$results$impact, imp$impact)
  # stronger omitted variable, larger impact
  expect_true(all(diff(imp$impact) > 0))

  out2 <- bootmakr(y ~ x + c, data = firms, treat = "x",
                   benchmark_covariates = "c", kd = c(0.5, 1, 2), reduce = FALSE,
                   reps = 20, seed = 1, progress = FALSE)
  imp2 <- out2$benchmark_strength$implied
  r_yz.x2 <- -imp2$r_yz.dx * sqrt((1 - r_yd^2) * (1 - imp2$r2dz.x)) + r_yd * imp2$r_dz.x
  expect_equal(imp2$impact, r_yz.x2 * imp2$r_dz.x, tolerance = 1e-12)
  expect_true(all(imp2$impact < imp$impact))
})

test_that("the printed block shows the impact column", {
  out <- bootmakr(y ~ x + c, data = firms, treat = "x",
                  benchmark_covariates = "c", reps = 20, seed = 1, progress = FALSE)
  txt <- capture.output(print(out))
  expect_true(any(grepl("Impact", txt)))
  expect_false(any(grepl("|r_dz.x|", txt, fixed = TRUE)))
})

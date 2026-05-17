# Shared fixtures for bootmakr tests.
# Skip every test in this suite if `sensemakr` is not installed (it provides
# the `darfur` dataset).
if (!requireNamespace("sensemakr", quietly = TRUE)) {
  skip_helper <- function() skip("sensemakr not installed")
} else {
  skip_helper <- function() invisible(NULL)
  utils::data("darfur", package = "sensemakr", envir = environment())
}

darfur_data <- if (exists("darfur")) darfur else NULL

darfur_formula <- peacefactor ~ directlyharmed + age + farmer_dar +
  herder_dar + pastvoted + hhsize_darfur + female + village

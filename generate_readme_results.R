## Regenerates README.md and the README figure in man/figures/ from README.Rmd.
## Since version 0.4.0 the code in README.Rmd runs when the file is rendered
## (a small simulated data set, a few seconds), so this is all that is needed.
## Run from the package root with the development version installed.
rmarkdown::render("README.Rmd",
                  output_format = rmarkdown::github_document(html_preview = FALSE))

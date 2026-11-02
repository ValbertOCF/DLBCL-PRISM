# Export the Shiny app to a static Shinylive site in docs/, preserving docs/data/.
# Run from the DLBCL-PRISM folder:  Rscript build.R
tmp <- tempfile("site")
shinylive::export("app", tmp)
keep <- c("data")
old <- setdiff(list.files("docs", all.files = FALSE), keep)
unlink(file.path("docs", old), recursive = TRUE)
file.copy(list.files(tmp, full.names = TRUE), "docs", recursive = TRUE)
file.create("docs/.nojekyll")
cat("Exported to docs/. Preview with: httpuv::runStaticServer('docs', port = 8008)\n")

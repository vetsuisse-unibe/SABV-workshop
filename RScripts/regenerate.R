#!/usr/bin/env Rscript
#' Regenerate the R-script versions of exercises 2a, 2b, 3a and 3b from their .qmd source.
#'
#' The .qmd files (in quarto_tutorials/) are the source of truth. The matching
#' 02a_/02b_/03a_/03b_*.R scripts are GENERATED copies: knitr::purl pulls the R code
#' out of the .qmd, keeping the prose as #' comments (documentation = 2), the same
#' style as RScripts/01_download_prepare_GTEx_liver.R.
#'
#' Re-run this whenever you edit any of those .qmd files, so the scripts stay in sync.
#'
#' Usage (from the repository root, or from inside RScripts/):
#'   Rscript RScripts/regenerate.R
#'
#' Note: exercise 01's script is maintained separately and is NOT regenerated here.

# Work whether this is launched from the repo root or from within RScripts/.
if (!dir.exists("quarto_tutorials") && dir.exists("../quarto_tutorials")) {
  setwd("..")
}
stopifnot(dir.exists("quarto_tutorials"), dir.exists("RScripts"))

tutorials <- c(
  "02a_workshop_simulation_exercises",
  "02b_workshop_gtex_liver_exploration",
  "03a_gtex_multitissue_interaction",
  "03b_gtex_multitissue_interaction_fullset"
)

for (name in tutorials) {
  knitr::purl(
    input         = file.path("quarto_tutorials", paste0(name, ".qmd")),
    output        = file.path("RScripts",         paste0(name, ".R")),
    documentation = 2,
    quiet         = TRUE
  )
  message("Regenerated RScripts/", name, ".R")
}

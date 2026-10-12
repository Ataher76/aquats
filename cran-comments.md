## Resubmission Summary
This is a resubmission of aquats version 1.2.2 addressing feedback from CRAN reviewer Leonore Hochhauser:

1. Acronym Expansion in DESCRIPTION:
   - All acronyms in the Description field have been expanded on first use (ANOVA, ANCOVA, MANOVA, PERMANOVA, SIMPER, NMDS, PCA, RDA, db-RDA, WQI).

2. Examples Formatting:
   - Removed \dontrun{} from R/summarize_obs.R.
   - Replaced example inputs with a fast, self-contained data frame so all examples execute deterministically in < 1 second.
   - Verified all 37 example files execute cleanly with devtools::run_examples(run_donttest = TRUE).

## Test environments
- local Windows 11, R 4.5.0

## R CMD check results
0 errors | 0 warnings | 0 notes

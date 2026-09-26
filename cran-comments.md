## Resubmission Summary (aquats 1.2.1)

This is a resubmission addressing comments from CRAN maintainers:

* Added formal methodology citations with persistent DOIs (`Authors (year) <doi:...>`) in the `Description` field of `DESCRIPTION`.
* Cleaned non-ASCII characters in `R/test_chisq.R`.
* Registered global variables (`stratum`) in `R/globals.R` to resolve `no visible binding` notes.
* Fixed empty link anchors in `README.md` and sanitized BibTeX URI formatting in vignette documentation.
* Added supervisor co-authorship metadata.

## Test environments
* local Windows 11, R 4.5.1
* win-builder (devel and release)

## R CMD check results
There were 0 errors | 0 warnings | 0 notes.

# aquats 1.2.1

## CRAN Maintenance & Resubmission Fixes

* **Methodology References**: Added primary literature citations formatted strictly as `Authors (year) <doi:...>` in the `Description` field of `DESCRIPTION` per CRAN maintainer request.
* **Non-ASCII Characters**: Removed non-ASCII encoding in `R/test_chisq.R` to ensure cross-platform compatibility.
* **Global Variables**: Declared `stratum` in `R/globals.R` to resolve the `no visible binding for global variable` check note in `plot_sankey()`.
* **Documentation & URLs**: Fixed empty badge links in `README.md` and sanitized BibTeX link syntax.
* **Authorship**: Updated metadata to include supervisor co-authorship.

---

# aquats 1.2.0

* Initial CRAN submission of `aquats`, providing publication-ready ecological, fisheries, and environmental data analysis and visualization workflows.

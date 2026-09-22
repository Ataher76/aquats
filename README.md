# aquats: Robust Ecological and Fisheries Data Analysis and Visualization

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)]()
[![Version](https://img.shields.io/badge/version-1.2.0-blue.svg)]()

**aquats** provides a streamlined, publication-grade analytical and visualization framework tailored for aquatic ecology, fisheries science, and environmental monitoring workflows.

---

## Key Modules & Capabilities

* **Parametric & Non-Linear Regressions:** Linear, allometric power ($y = a \cdot x^b$), exponential, logarithmic, polynomial, and cubic curves with auto-generated algebraic labels, $R^2$, and $p$-values.
* **Univariate & Multivariate Statistics:** Streamlined wrappers for one-way/factorial ANOVA, ANCOVA, MANOVA, and Generalized Linear Models (GLMs) for count and continuous data.
* **Multivariate Ordination & Community Ecology:** Publication-ready NMDS, PCA, RDA, distance-based RDA (db-RDA), hierarchical cluster dendrograms, and integrated Mantel test network heatmaps.
* **Ecological Diversity & Accumulation:** Alpha diversity profiling (Shannon, Simpson, Chao1, Pielou evenness) and species accumulation curves.
* **Environmental & Ecotoxicological Assessment:** Weighted Arithmetic Water Quality Index (WAWQI) and comprehensive heavy metal contamination metrics (CF, CD, PLI, HPI, HEI).
* **Longitudinal Time-Series Tracking:** Multi-station temporal trends with LOESS/linear fits and 95% confidence ribbons.

---

## Installation

You can install the development version of `aquats` directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Ataher76/aquats")

```

---

## Quick Start

### 1. Non-Linear Allometric Regression

```r
library(aquats)

# Fit an allometric power curve with 95% confidence intervals
res <- regression_analysis(
  data = hilsa_master,
  x_var = "TL_cm",
  y_var = "BW_g",
  fit_type = "power",
  interval = "confidence",
  point_size = 2.5
)

# Display publication-ready plot
print(res$Plot)

```

### 2. Community Alpha Diversity Profiling

```r
# Calculate alpha diversity indices using bundled community data
div_res <- calc_diversity(
  data = community_master
)

head(div_res)

```

## Author & Maintainer

* **Ataher Ali** – Maintainer & Author
* ORCID: [0009-0008-8780-0903](https://orcid.org/0009-0008-8780-0903)
* Email: `ataher.cu.ms@gmail.com`

## License

This package is licensed under the **MIT License**.


#' Comprehensive Master Dataset for Hilsa Shad Ecology and Statistics
#'
#' A comprehensive, unified dataset containing morphological measurements, environmental
#' parameters, catch records, and multi-factor experimental groupings for Hilsa shad
#' (Tenualosa ilisha) across various habitats, seasons, and size classes. Designed to
#' power all parametric, non-parametric, regression, and multivariate functions in the package.
#'
#' @format A data frame with 120 rows and 14 variables:
#' \describe{
#'   \item{Habitat}{Factor; sampling environment (Marine, Estuary, River).}
#'   \item{Season}{Factor; fishing season (Monsoon, Dry).}
#'   \item{Size_Class}{Factor; growth stage (Juvenile, Adult).}
#'   \item{Weight_g}{Numeric; body weight in grams.}
#'   \item{Total_Length_cm}{Numeric; total length in centimeters.}
#'   \item{Body_Depth_cm}{Numeric; maximum body depth in centimeters.}
#'   \item{Head_Length_cm}{Numeric; head length in centimeters.}
#'   \item{Fin_Length_cm}{Numeric; pectoral fin length in centimeters.}
#'   \item{Fishing_Hours}{Numeric; fishing effort in hours.}
#'   \item{Catch_Count}{Integer; number of fish caught (count response).}
#'   \item{SST_C}{Numeric; Sea Surface Temperature in degrees Celsius.}
#'   \item{Salinity_ppt}{Numeric; water salinity in parts per thousand.}
#'   \item{Depth_m}{Numeric; water depth in meters.}
#'   \item{Catch_kg}{Numeric; total fish catch in kilograms.}
#' }
"hilsa_master"


#' Multi-Trophic Community Abundance Master List
#'
#' A curated list containing species abundance matrices for phytoplankton, zooplankton,
#' and benthic macroinvertebrates, designed for community ecology analyses, PERMANOVA,
#' SIMPER, and multi-matrix Mantel network visualizations.
#'
#' @format A named list with 3 data frames:
#' \describe{
#'   \item{Phytoplankton}{Data frame of diatom, cyanobacteria, and green algae counts.}
#'   \item{Zooplankton}{Data frame of copepod, cladoceran, and rotifer counts.}
#'   \item{Benthos}{Data frame of chironomid, oligochaete, and mollusc counts.}
#' }
"community_master"




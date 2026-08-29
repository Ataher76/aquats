#' Body Weights of Hilsa Shad Across Habitats
#'
#' A simulated dataset containing the body weights of Hilsa shad
#' (Tenualosa ilisha) sampled from three distinct ecological environments:
#' Marine, Estuary, and River. This dataset is designed to demonstrate
#' one-way ANOVA and variance visualization.
#'
#' @format A data frame with 90 rows and 2 variables:
#' \describe{
#'   \item{Habitat}{A factor representing the sampling environment (Marine, Estuary, River).}
#'   \item{Weight_g}{A numeric vector representing the body weight of the fish in grams.}
#' }
"hilsa_weight"


#' Multi-Factor Body Weights of Hilsa Shad
#'
#' A simulated dataset containing the body weights of Hilsa shad
#' (Tenualosa ilisha) across different habitats and fishing seasons
#' to demonstrate two-way ANOVA interaction effects.
#'
#' @format A data frame with 90 rows and 3 variables:
#' \describe{
#'   \item{Habitat}{A factor representing the environment (Marine, Estuary, River).}
#'   \item{Season}{A factor representing the season (Monsoon, Dry).}
#'   \item{Weight_g}{A numeric vector representing fish body weight in grams.}
#' }
"hilsa_two_way"



#' Comprehensive Body Weights of Hilsa Shad (Three-Way)
#'
#' A simulated dataset containing body weights of Hilsa shad (Tenualosa ilisha)
#' across multiple habitats, fishing seasons, and size classes to demonstrate
#' three-way ANOVA interaction models and faceted visualization.
#'
#' @format A data frame with 120 rows and 4 variables:
#' \describe{
#'   \item{Habitat}{A factor representing the environment (Marine, Estuary, River).}
#'   \item{Season}{A factor representing the season (Monsoon, Dry).}
#'   \item{Size_Class}{A factor representing growth stage (Juvenile, Adult).}
#'   \item{Weight_g}{A numeric vector representing fish body weight in grams.}
#' }
"hilsa_three_way"


#' Morphometric Measurements of Hilsa Shad
#'
#' A simulated dataset containing multiple continuous morphological traits
#' of Hilsa shad (Tenualosa ilisha) across different aquatic habitats to
#' demonstrate Multivariate Analysis of Variance (MANOVA).
#'
#' @format A data frame with 90 rows and 4 variables:
#' \describe{
#'   \item{Habitat}{A factor representing the environment (Marine, Estuary, River).}
#'   \item{Body_Depth_cm}{A numeric vector representing maximum body depth in centimeters.}
#'   \item{Head_Length_cm}{A numeric vector representing head length in centimeters.}
#'   \item{Fin_Length_cm}{A numeric vector representing pectoral fin length in centimeters.}
#' }
"hilsa_morphology"


#' Length-Weight ANCOVA Dataset for Hilsa Shad
#'
#' A simulated dataset containing body weight, total length, and habitat types
#' of Hilsa shad (Tenualosa ilisha) to demonstrate Analysis of Covariance (ANCOVA).
#'
#' @format A data frame with 90 rows and 3 variables:
#' \describe{
#'   \item{Habitat}{A factor representing the environment (Marine, Estuary, River).}
#'   \item{Total_Length_cm}{A numeric vector representing total fish length in centimeters (covariate).}
#'   \item{Weight_g}{A numeric vector representing fish body weight in grams (response).}
#' }
"hilsa_ancova"



#' Morphometric Regression Dataset for Hilsa Shad
#'
#' A simulated dataset containing length and weight measurements of Hilsa shad
#' (Tenualosa ilisha) across various habitats for evaluating linear, logarithmic,
#' and polynomial regression models.
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{Total_Length_cm}{A numeric vector representing total fish length in centimeters.}
#'   \item{Weight_g}{A numeric vector representing fish body weight in grams.}
#'   \item{Habitat}{A factor representing the environment (River, Estuary, Marine).}
#' }
"hilsa_regression"


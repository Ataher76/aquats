#' Analysis of Covariance (ANCOVA) with Adjusted Means Visualization
#'
#' This function performs an Analysis of Covariance to test group differences in a
#' continuous response variable while controlling for a continuous covariate. It extracts
#' the ANCOVA table, estimated marginal means (EMMs), and generates a publication-ready
#' regression scatter plot with parallel trend lines.
#'
#' @param data A data frame in long format.
#' @param response_var Character; the name of the continuous numeric response variable.
#' @param factor_var Character; the name of the categorical independent variable (groups).
#' @param covariate_var Character; the name of the continuous covariate variable.
#' @param factor_levels Optional character vector; custom order for the factor levels.
#' @param color_palette Character; name of a color palette from RColorBrewer. Default is "Set1".
#'
#' @return A list containing the ANCOVA model summary, EMMs, and the ggplot object.
#'
#' @import ggplot2 emmeans RColorBrewer dplyr
#' @importFrom stats lm anova as.formula
#' @importFrom dplyr %>%
#' @export
ancova_analysis <- function(data,
                            response_var,
                            factor_var,
                            covariate_var,
                            factor_levels = NULL,
                            color_palette = "Set1") {

  # 1. Input Validation & Factor Ordering
  required_cols <- c(factor_var, response_var, covariate_var)
  if (!all(required_cols %in% names(data))) stop("Specified columns not found in data.")
  if (!is.numeric(data[[response_var]])) stop(paste(response_var, "must be numeric."))
  if (!is.numeric(data[[covariate_var]])) stop(paste(covariate_var, "must be numeric."))

  if (!is.null(factor_levels)) {
    data[[factor_var]] <- factor(data[[factor_var]], levels = factor_levels)
  } else {
    data[[factor_var]] <- as.factor(data[[factor_var]])
  }

  # 2. Fit ANCOVA Model (Additive model assuming homogeneity of slopes)
  formula_str <- as.formula(paste(response_var, "~", covariate_var, "+", factor_var))
  ancova_model <- stats::lm(formula_str, data = data)
  ancova_table <- stats::anova(ancova_model)

  # 3. Calculate Estimated Marginal Means (EMMs) adjusted for covariate
  emms <- emmeans::emmeans(ancova_model, specs = as.formula(paste("~", factor_var)))

  # 4. Visualization: Scatter plot with regression trend lines
  p <- ggplot(data, aes(x = .data[[covariate_var]], y = .data[[response_var]], color = .data[[factor_var]])) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
    theme_bw() +
    labs(
      x = covariate_var,
      y = response_var,
      color = factor_var,
      fill = factor_var,
      title = "Analysis of Covariance (ANCOVA)"
    ) +
    scale_color_brewer(palette = color_palette) +
    scale_fill_brewer(palette = color_palette) +
    theme(plot.title = element_text(hjust = 0.5))

  return(list(
    ANCOVA_Table = ancova_table,
    Adjusted_Means = emms,
    Plot = p
  ))
}


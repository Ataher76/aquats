#' Analysis of Covariance (ANCOVA) with Adjusted Means Visualization
#'
#' This function performs an Analysis of Covariance to test group differences in a
#' continuous response variable while controlling for a continuous covariate. It extracts
#' the ANCOVA table, estimated marginal means (EMMs), and generates a publication-ready
#' regression scatter plot with parallel trend lines matching the additive model.
#'
#' @param data A data frame in long format.
#' @param response_var Character; the name of the continuous numeric response variable.
#' @param factor_var Character; the name of the categorical independent variable (groups).
#' @param covariate_var Character; the name of the continuous covariate variable.
#' @param factor_levels Optional character vector; custom order for the factor levels.
#' @param facet Logical. If \code{TRUE}, panels are split by factor levels using \code{facet_wrap} (default: \code{FALSE}).
#' @param color_palette Character; name of a color palette from RColorBrewer. Default is "Set1".
#'
#' @return A list containing the ANCOVA model summary, EMMs, and the ggplot object.
#'
#' @export
#' @import ggplot2 emmeans RColorBrewer
#' @importFrom stats lm anova as.formula predict
#' @importFrom rlang .data
ancova_analysis <- function(data,
                            response_var,
                            factor_var,
                            covariate_var,
                            factor_levels = NULL,
                            facet = FALSE,
                            color_palette = "Set1") {

  # Standard Defensive Ingestion Guard
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  # 1. Input Validation & Defensive NA Cleaning
  required_cols <- c(factor_var, response_var, covariate_var)
  if (!all(required_cols %in% names(data))) {
    stop("Specified response, factor, or covariate columns not found in data.")
  }
  if (!is.numeric(data[[response_var]])) stop(paste(response_var, "must be numeric."))
  if (!is.numeric(data[[covariate_var]])) stop(paste(covariate_var, "must be numeric."))

  valid_rows <- !is.na(data[[factor_var]]) & !is.na(data[[response_var]]) & !is.na(data[[covariate_var]])
  clean_data <- data[valid_rows, , drop = FALSE]
  na_dropped <- nrow(data) - nrow(clean_data)

  if (na_dropped > 0) {
    message(sprintf("Note: Automatically removed %d row(s) containing NA values in target variables.", na_dropped))
  }

  if (nrow(clean_data) < 5) {
    stop("ANCOVA requires at least 5 complete observations.")
  }

  if (!is.null(factor_levels)) {
    clean_data[[factor_var]] <- factor(clean_data[[factor_var]], levels = factor_levels)
  } else {
    clean_data[[factor_var]] <- as.factor(clean_data[[factor_var]])
  }

  if (nlevels(clean_data[[factor_var]]) < 2) {
    stop("The grouping factor must contain at least 2 distinct levels for ANCOVA.")
  }

  # 2. Fit ANCOVA Model (Additive model assuming homogeneity of slopes)
  formula_str <- as.formula(paste(response_var, "~", covariate_var, "+", factor_var))
  ancova_model <- stats::lm(formula_str, data = clean_data)
  ancova_table <- stats::anova(ancova_model)

  clean_data$.fitted <- stats::predict(ancova_model)

  # 3. Calculate Estimated Marginal Means (EMMs) adjusted for covariate
  emms <- emmeans::emmeans(ancova_model, specs = as.formula(paste("~", factor_var)))

  # 4. Visualization
  p <- ggplot2::ggplot(
    clean_data,
    ggplot2::aes(x = .data[[covariate_var]], y = .data[[response_var]], color = .data[[factor_var]])
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 2.2) +
    ggplot2::geom_line(ggplot2::aes(y = .data$.fitted), linewidth = 1.1) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = covariate_var,
      y = response_var,
      color = factor_var,
      fill = factor_var,
      title = "Analysis of Covariance (ANCOVA)"
    ) +
    ggplot2::scale_color_brewer(palette = color_palette) +
    ggplot2::scale_fill_brewer(palette = color_palette) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      legend.position = if (facet) "none" else "right",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    )

  # Optional Facet Wrap
  if (facet) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", factor_var)), scales = "fixed")
  }

  return(list(
    ANCOVA_Table = ancova_table,
    Adjusted_Means = emms,
    Plot = p
  ))
}

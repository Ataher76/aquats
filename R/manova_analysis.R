#' Multivariate Analysis of Variance (MANOVA) with Alternative Visualizations
#'
#' This function performs a one-way MANOVA across multiple continuous response
#' variables grouped by a categorical factor, extracts multivariate test statistics,
#' univariate ANOVA breakdowns, and generates either faceted boxplots or a Canonical
#' Discriminant Analysis (LDA) multivariate scatter plot.
#'
#' @param data A data frame in long format.
#' @param response_vars Character vector; the names of the continuous numeric response variables.
#' @param factor_var Character; the name of the categorical independent variable.
#' @param factor_levels Optional character vector; custom order for the factor levels.
#' @param plot_type Character; type of visualization: "boxplot" or "lda". Default is "boxplot".
#' @param color_palette Character; name of a color palette from RColorBrewer. Default is "Set1".
#'
#' @return A list containing the MANOVA test summary, univariate ANOVA summaries, and the ggplot object.
#'
#' @import ggplot2 RColorBrewer dplyr tidyr MASS
#' @importFrom stats manova lm as.formula summary.aov predict
#' @importFrom dplyr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom MASS lda
#' @export
manova_analysis <- function(data,
                            response_vars,
                            factor_var,
                            factor_levels = NULL,
                            plot_type = "boxplot",
                            color_palette = "Set1") {

  # 1. Input Validation & Factor Ordering
  required_cols <- c(factor_var, response_vars)
  if (!all(required_cols %in% names(data))) stop("Specified columns not found in data.")
  for (v in response_vars) {
    if (!is.numeric(data[[v]])) stop(paste(v, "must be a numeric variable."))
  }

  if (!is.null(factor_levels)) {
    data[[factor_var]] <- factor(data[[factor_var]], levels = factor_levels)
  } else {
    data[[factor_var]] <- as.factor(data[[factor_var]])
  }

  # 2. Fit Multivariate Linear Model & MANOVA
  lhs_formula <- paste("cbind(", paste(response_vars, collapse = ", "), ")")
  formula_str <- as.formula(paste(lhs_formula, "~", factor_var))

  mlm_model <- stats::lm(formula_str, data = data)
  manova_res <- stats::manova(mlm_model)
  manova_summary <- summary(manova_res, test = "Pillai")
  univariate_anovas <- summary.aov(manova_res)

  # 3. Plot Routing
  if (plot_type == "boxplot") {
    # Faceted Boxplot Visualization
    plot_data <- data %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(response_vars),
        names_to = "Morphometric_Trait",
        values_to = "Measurement"
      )

    p <- ggplot(plot_data, aes(x = .data[[factor_var]], y = Measurement, fill = .data[[factor_var]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_point(position = position_jitter(width = 0.1), size = 1.2, alpha = 0.4, show.legend = FALSE) +
      facet_wrap(~ Morphometric_Trait, scales = "free_y") +
      theme_bw() +
      labs(
        x = factor_var,
        y = "Measurement Value",
        fill = factor_var,
        title = "Multivariate Analysis of Variance (MANOVA) - Traits"
      ) +
      scale_fill_brewer(palette = color_palette) +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill = "lightgray"))

  } else if (plot_type == "lda") {
    # Canonical Discriminant Analysis (LDA) Scatter Plot
    lda_formula <- as.formula(paste(factor_var, "~", paste(response_vars, collapse = " + ")))
    lda_model <- MASS::lda(lda_formula, data = data)
    lda_preds <- predict(lda_model)

    lda_df <- data.frame(
      Group = data[[factor_var]],
      LD1 = lda_preds$x[, 1],
      LD2 = if (ncol(lda_preds$x) > 1) lda_preds$x[, 2] else rep(0, nrow(data))
    )

    # Determine axes labels based on available discriminant dimensions
    x_lab <- paste0("Canonical Dimension 1 (LD1)")
    y_lab <- if (ncol(lda_preds$x) > 1) "Canonical Dimension 2 (LD2)" else "Dimension 2 (N/A)"

    p <- ggplot(lda_df, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
      geom_point(size = 2.5, alpha = 0.8) +
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE) +
      theme_bw() +
      labs(
        x = x_lab,
        y = y_lab,
        color = factor_var,
        fill = factor_var,
        title = "Canonical Discriminant Analysis (MANOVA Space)"
      ) +
      scale_color_brewer(palette = color_palette) +
      scale_fill_brewer(palette = color_palette) +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    stop("Invalid plot_type. Choose either 'boxplot' or 'lda'.")
  }

  return(list(
    MANOVA_Summary = manova_summary,
    Univariate_ANOVAs = univariate_anovas,
    Plot = p
  ))
}

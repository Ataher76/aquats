#' Multivariate Analysis of Variance (MANOVA) with Alternative Visualizations
#'
#' Performs a one-way MANOVA across multiple continuous response variables grouped
#' by a categorical factor, extracts multivariate test statistics (Pillai's Trace),
#' provides univariate ANOVA breakdowns, and generates either faceted trait boxplots
#' or a Canonical Discriminant Analysis (LDA) multivariate ordination scatter plot.
#'
#' @param data A data frame or matrix in long format.
#' @param response_vars Character vector; the names of the continuous numeric response variables.
#' @param factor_var Character; the name of the categorical independent variable.
#' @param factor_levels Optional character vector; custom order for the factor levels.
#' @param plot_type Character; type of visualization: "boxplot" or "lda". Default is "boxplot".
#' @param color_palette Character; name of a color palette from RColorBrewer. Default is "Set1".
#'
#' @return A list containing:
#'   \item{MANOVA_Summary}{Multivariate test summary based on Pillai's trace.}
#'   \item{Univariate_ANOVAs}{Univariate ANOVA tables for each response variable.}
#'   \item{Plot}{Publication-ready ggplot object.}
#' @export
#'
#' @import ggplot2 RColorBrewer
#' @importFrom stats manova lm as.formula summary.aov predict complete.cases
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr all_of
#' @importFrom MASS lda
#' @importFrom rlang .data
manova_analysis <- function(data,
                            response_vars,
                            factor_var,
                            factor_levels = NULL,
                            plot_type = "boxplot",
                            color_palette = "Set1") {

  # Standard Defensive Ingestion Guard
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  # 1. Input Validation & Defensive NA Cleaning
  required_cols <- c(factor_var, response_vars)
  if (!all(required_cols %in% names(data))) {
    missing_cols <- required_cols[!required_cols %in% names(data)]
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }

  for (v in response_vars) {
    if (!is.numeric(data[[v]])) stop(paste(v, "must be a numeric variable."))
  }

  # Isolate valid rows exclusively for target variables
  clean_subset <- data[, required_cols, drop = FALSE]
  valid_rows <- stats::complete.cases(clean_subset)
  clean_data <- data[valid_rows, , drop = FALSE]

  na_dropped <- sum(!valid_rows)
  if (na_dropped > 0) {
    message(sprintf("Note: Automatically removed %d row(s) containing NA values in target MANOVA variables.", na_dropped))
  }

  if (nrow(clean_data) < (length(response_vars) + 2)) {
    stop("Insufficient complete observations to estimate MANOVA covariance matrices.")
  }

  # Factor Ordering / Coercion
  if (!is.null(factor_levels)) {
    clean_data[[factor_var]] <- factor(clean_data[[factor_var]], levels = factor_levels)
  } else {
    clean_data[[factor_var]] <- as.factor(clean_data[[factor_var]])
  }

  if (nlevels(clean_data[[factor_var]]) < 2) {
    stop("The grouping factor must contain at least 2 distinct levels for MANOVA.")
  }

  # 2. Fit Multivariate Linear Model & MANOVA
  lhs_formula <- paste0("cbind(", paste(response_vars, collapse = ", "), ")")
  formula_str <- stats::as.formula(paste(lhs_formula, "~", factor_var))

  mlm_model <- stats::lm(formula_str, data = clean_data)
  manova_res <- stats::manova(mlm_model)
  manova_summary <- summary(manova_res, test = "Pillai")
  univariate_anovas <- stats::summary.aov(manova_res)

  # 3. Plot Routing
  if (plot_type == "boxplot") {
    # Faceted Boxplot Visualization across all continuous response traits
    plot_data <- tidyr::pivot_longer(
      clean_data,
      cols = dplyr::all_of(response_vars),
      names_to = "Response_Trait",
      values_to = "Measurement"
    )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data[[factor_var]], y = .data$Measurement, fill = .data[[factor_var]])
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = 0.15),
        size = 1.2, alpha = 0.45, color = "black", show.legend = FALSE
      ) +
      ggplot2::facet_wrap(~ Response_Trait, scales = "free_y") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::scale_fill_brewer(palette = color_palette) +
      ggplot2::labs(
        x = factor_var,
        y = "Measurement Value",
        fill = factor_var,
        title = "Multivariate Analysis of Variance (MANOVA) - Trait Breakdown"
      ) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black", face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
      )

  } else if (plot_type == "lda") {
    # Canonical Discriminant Analysis (LDA) Ordination Plot
    lda_formula <- stats::as.formula(paste(factor_var, "~", paste(response_vars, collapse = " + ")))
    lda_model <- MASS::lda(lda_formula, data = clean_data)
    lda_preds <- stats::predict(lda_model)

    num_dims <- ncol(lda_preds$x)

    if (num_dims >= 2) {
      lda_df <- data.frame(
        Group = clean_data[[factor_var]],
        LD1 = lda_preds$x[, 1],
        LD2 = lda_preds$x[, 2]
      )

      p <- ggplot2::ggplot(lda_df, ggplot2::aes(x = .data$LD1, y = .data$LD2, color = .data$Group, fill = .data$Group)) +
        ggplot2::geom_point(size = 2.8, alpha = 0.8) +
        ggplot2::stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(
          x = "Canonical Dimension 1 (LD1)",
          y = "Canonical Dimension 2 (LD2)",
          color = factor_var,
          fill = factor_var,
          title = "Canonical Discriminant Analysis (MANOVA Space)"
        ) +
        ggplot2::scale_color_brewer(palette = color_palette) +
        ggplot2::scale_fill_brewer(palette = color_palette) +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
          axis.title = ggplot2::element_text(face = "bold"),
          axis.text = ggplot2::element_text(color = "black", face = "bold"),
          legend.title = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
        )
    } else {
      # Fallback for 2-group factor: 1D Density Ridge / Histogram along LD1
      lda_df <- data.frame(
        Group = clean_data[[factor_var]],
        LD1 = lda_preds$x[, 1]
      )

      p <- ggplot2::ggplot(lda_df, ggplot2::aes(x = .data$LD1, fill = .data$Group, color = .data$Group)) +
        ggplot2::geom_density(alpha = 0.4, linewidth = 1) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(
          x = "Canonical Discriminant Axis 1 (LD1)",
          y = "Density",
          fill = factor_var,
          color = factor_var,
          title = "Canonical Discriminant Separation (2 Groups)"
        ) +
        ggplot2::scale_color_brewer(palette = color_palette) +
        ggplot2::scale_fill_brewer(palette = color_palette) +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
          axis.title = ggplot2::element_text(face = "bold"),
          axis.text = ggplot2::element_text(color = "black", face = "bold"),
          legend.title = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
        )
    }

  } else {
    stop("Invalid plot_type. Choose either 'boxplot' or 'lda'.")
  }

  return(list(
    MANOVA_Summary = manova_summary,
    Univariate_ANOVAs = univariate_anovas,
    Plot = p
  ))
}

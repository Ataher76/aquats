#' Three-Way Analysis of Variance (ANOVA) with Faceted Interaction Visualization
#'
#' Performs a three-way factorial ANOVA with full interaction effects, calculates
#' Tukey's HSD post-hoc tests, extracts Estimated Marginal Means (EMMs), and generates
#' faceted publication-ready interaction plots with compact letter displays.
#' Automatically handles NA values safely, converts matrix inputs, and coerces categorical variables.
#'
#' @param data A data frame or matrix in long format.
#' @param factor1_var Character; the first categorical independent variable (x-axis).
#' @param factor2_var Character; the second categorical independent variable (groups/fill).
#' @param factor3_var Character; the third categorical independent variable (facets).
#' @param numeric_var Character; the continuous response variable.
#' @param factor1_levels Optional character vector; custom order for the first factor.
#' @param factor2_levels Optional character vector; custom order for the second factor.
#' @param factor3_levels Optional character vector; custom order for the third factor.
#' @param plot_type Character; "boxplot", "barplot", or "pointrange". Default is "boxplot".
#' @param error_type Character; error bar type: "se" (Standard Error) or "sd" (Standard Deviation).
#' @param y_limits Optional numeric vector of length 2; explicitly sets the Y-axis limits.
#' @param add_jitter Logical; if TRUE, adds jittered points to boxplots. Default is TRUE.
#' @param show_mean Logical; if TRUE, displays mean points inside boxplots. Default is TRUE.
#' @param mean_color Character; color for mean markers. Default is "darkred".
#' @param color_palette Character; valid RColorBrewer palette name. Default is "Set1".
#'
#' @return A list containing:
#'   \item{ANOVA_Summary}{Global ANOVA summary table with 3-way interaction terms.}
#'   \item{TukeyHSD}{Tukey HSD pairwise comparison results.}
#'   \item{EMMs}{Estimated Marginal Means object.}
#'   \item{Summary_Table}{Aggregated means, SDs, SEs, and confidence limits.}
#'   \item{Plot}{The publication-ready faceted ggplot object.}
#' @export
#'
#' @import ggplot2 emmeans multcomp RColorBrewer
#' @importFrom stats aov as.formula sd TukeyHSD complete.cases
#' @importFrom dplyr group_by summarise mutate across all_of
#' @importFrom rlang .data
three_way_anova <- function(data,
                            factor1_var,
                            factor2_var,
                            factor3_var,
                            numeric_var,
                            factor1_levels = NULL,
                            factor2_levels = NULL,
                            factor3_levels = NULL,
                            plot_type = "boxplot",
                            error_type = "se",
                            y_limits = NULL,
                            add_jitter = TRUE,
                            show_mean = TRUE,
                            mean_color = "darkred",
                            color_palette = "Set1") {

  # 1. Defensive Ingestion
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  required_cols <- c(factor1_var, factor2_var, factor3_var, numeric_var)
  if (!all(required_cols %in% names(data))) {
    missing_cols <- required_cols[!required_cols %in% names(data)]
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(data[[numeric_var]])) {
    stop(paste(numeric_var, "must be numeric."))
  }

  # Complete-case filtering
  valid_rows <- stats::complete.cases(data[, required_cols, drop = FALSE])
  clean_data <- data[valid_rows, , drop = FALSE]
  na_dropped <- nrow(data) - nrow(clean_data)
  if (na_dropped > 0) {
    message(sprintf("Note: Automatically removed %d row(s) containing NA values in target variables.", na_dropped))
  }

  if (nrow(clean_data) < 8) {
    stop("Three-Way ANOVA requires at least 8 valid observations.")
  }

  # Factor Ordering & Coercion
  if (!is.null(factor1_levels)) {
    clean_data[[factor1_var]] <- factor(clean_data[[factor1_var]], levels = factor1_levels)
  } else {
    clean_data[[factor1_var]] <- as.factor(clean_data[[factor1_var]])
  }

  if (!is.null(factor2_levels)) {
    clean_data[[factor2_var]] <- factor(clean_data[[factor2_var]], levels = factor2_levels)
  } else {
    clean_data[[factor2_var]] <- as.factor(clean_data[[factor2_var]])
  }

  if (!is.null(factor3_levels)) {
    clean_data[[factor3_var]] <- factor(clean_data[[factor3_var]], levels = factor3_levels)
  } else {
    clean_data[[factor3_var]] <- as.factor(clean_data[[factor3_var]])
  }

  # 2. Factorial Model Fitting
  formula_str <- stats::as.formula(paste(numeric_var, "~", factor1_var, "*", factor2_var, "*", factor3_var))
  aov_model <- stats::aov(formula_str, data = clean_data)
  aov_summary <- summary(aov_model)
  tukey_res <- stats::TukeyHSD(aov_model)

  # 3. Estimated Marginal Means & CLD
  emmean <- emmeans::emmeans(
    aov_model,
    specs = stats::as.formula(paste("~", factor1_var, "*", factor2_var, "*", factor3_var))
  )
  emmean_cld <- multcomp::cld(emmean, Letters = letters)
  emmean_cld$.group <- trimws(emmean_cld$.group)

  # 4. Summary Statistics for Error Bars & Positioning
  sum_data <- clean_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(factor1_var, factor2_var, factor3_var)))) %>%
    dplyr::summarise(
      mean_val = mean(.data[[numeric_var]], na.rm = TRUE),
      sd_val = stats::sd(.data[[numeric_var]], na.rm = TRUE),
      se_val = sd_val / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      error_margin = if (error_type == "sd") sd_val else se_val,
      upper_limit = mean_val + error_margin
    )

  emmean_cld <- merge(emmean_cld, sum_data, by = c(factor1_var, factor2_var, factor3_var))

  # 5. Base Plot Initialization with Faceting
  p <- ggplot2::ggplot(
    clean_data,
    ggplot2::aes(x = .data[[factor1_var]], y = .data[[numeric_var]], fill = .data[[factor2_var]])
  ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", factor3_var))) +
    ggplot2::scale_fill_brewer(palette = color_palette) +
    ggplot2::scale_color_brewer(palette = color_palette) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    ) +
    ggplot2::labs(
      x = factor1_var,
      y = numeric_var,
      fill = factor2_var,
      color = factor2_var,
      title = "Three-Way ANOVA Interaction"
    )

  # 6. Plot Type Geometries
  if (plot_type == "boxplot") {
    p <- p +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75, position = ggplot2::position_dodge(0.8))

    if (add_jitter) {
      p <- p + ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
        color = "black", size = 1.0, alpha = 0.4, show.legend = FALSE
      )
    }

    if (show_mean) {
      p <- p + ggplot2::stat_summary(
        fun = mean, geom = "point", shape = 18, size = 3,
        color = mean_color, position = ggplot2::position_dodge(0.8), show.legend = FALSE
      )
    }

  } else if (plot_type == "barplot") {
    p <- ggplot2::ggplot(
      sum_data,
      ggplot2::aes(x = .data[[factor1_var]], y = .data$mean_val, fill = .data[[factor2_var]])
    ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", factor3_var))) +
      ggplot2::geom_col(alpha = 0.8, color = "black", position = ggplot2::position_dodge(0.8)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$mean_val - .data$error_margin, ymax = .data$upper_limit),
        width = 0.2, position = ggplot2::position_dodge(0.8)
      ) +
      ggplot2::scale_fill_brewer(palette = color_palette) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black", face = "bold"),
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
      ) +
      ggplot2::labs(
        x = factor1_var,
        y = numeric_var,
        fill = factor2_var,
        title = paste("Three-Way ANOVA (Mean \u00B1", toupper(error_type), ")")
      )

  } else if (plot_type == "pointrange") {
    p <- ggplot2::ggplot(
      sum_data,
      ggplot2::aes(x = .data[[factor1_var]], y = .data$mean_val, color = .data[[factor2_var]])
    ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", factor3_var))) +
      ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(0.4)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$mean_val - .data$error_margin, ymax = .data$upper_limit),
        width = 0.2, position = ggplot2::position_dodge(0.4)
      ) +
      ggplot2::scale_color_brewer(palette = color_palette) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black", face = "bold"),
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
      ) +
      ggplot2::labs(
        x = factor1_var,
        y = numeric_var,
        color = factor2_var,
        title = paste("Three-Way ANOVA (Mean \u00B1", toupper(error_type), ")")
      )
  }

  # 7. Post-Hoc Letter Display Placement
  y_nudge <- max(clean_data[[numeric_var]], na.rm = TRUE) * 0.05
  if (plot_type == "boxplot") {
    max_y <- max(clean_data[[numeric_var]], na.rm = TRUE)
    p <- p + ggplot2::geom_text(
      data = emmean_cld,
      ggplot2::aes(x = .data[[factor1_var]], y = max_y + y_nudge, label = .data$.group, group = .data[[factor2_var]]),
      position = ggplot2::position_dodge(0.8), inherit.aes = FALSE, size = 3.6, fontface = "bold"
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = emmean_cld,
      ggplot2::aes(x = .data[[factor1_var]], y = .data$upper_limit + y_nudge, label = .data$.group, group = .data[[factor2_var]]),
      position = ggplot2::position_dodge(0.8), inherit.aes = FALSE, size = 3.6, fontface = "bold"
    )
  }

  # 8. Custom Y Limits
  if (!is.null(y_limits)) {
    if (length(y_limits) == 2 && is.numeric(y_limits)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    } else {
      warning("y_limits must be a numeric vector of length 2. Ignoring custom limits.")
    }
  }

  return(list(
    ANOVA_Summary = aov_summary,
    TukeyHSD = tukey_res,
    EMMs = emmean,
    Summary_Table = sum_data,
    Plot = p
  ))
}

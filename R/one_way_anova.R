#' One-Way Analysis of Variance (ANOVA) with Visualization
#'
#' This function performs a one-way ANOVA, calculates Tukey's HSD post-hoc test,
#' extracts Estimated Marginal Means (EMMs), and generates a publication-ready plot.
#'
#' @param data A data frame in long format.
#' @param factor_var Character; the name of the categorical independent variable.
#' @param numeric_var Character; the name of the continuous response variable.
#' @param factor_levels Optional character vector; specifies the exact order of the categorical levels on the x-axis.
#' @param plot_type Character; type of plot: "boxplot", "barplot", or "pointrange". Default is "boxplot".
#' @param error_type Character; error bar type for bar/pointrange plots: "se" (Standard Error) or "sd" (Standard Deviation).
#' @param sig_display Character; display significance via "letters" (compact letter display) or "stars" (p-value brackets).
#' @param y_limits Optional numeric vector of length 2; explicitly sets the Y-axis limits (e.g., c(0, 1500)).
#' @param add_jitter Logical; if TRUE, adds jittered points to the boxplot.
#' @param show_mean Logical; if TRUE, displays a mean point inside the boxplot.
#' @param mean_color Character; color for the mean point in the boxplot.
#' @param color_palette Character; name of a color palette from RColorBrewer.
#'
#' @return A list containing the ANOVA summary, Tukey HSD results, EMMs, a clean summary table, and the ggplot object.
#'
#' @import ggplot2 emmeans multcomp RColorBrewer ggpubr dplyr
#' @importFrom stats aov as.formula sd
#' @export
one_way_anova <- function(data,
                          factor_var,
                          numeric_var,
                          factor_levels = NULL,
                          plot_type = "boxplot",
                          error_type = "se",
                          sig_display = "letters",
                          y_limits = NULL,
                          add_jitter = TRUE,
                          show_mean = TRUE,
                          mean_color = "darkred",
                          color_palette = "Set1") {

  # 1. Input Validation & Factor Ordering
  if (!factor_var %in% names(data) || !numeric_var %in% names(data)) stop("Specified columns not found in data.")
  if (!is.numeric(data[[numeric_var]])) stop(paste(numeric_var, "must be numeric."))

  if (!is.null(factor_levels)) {
    data[[factor_var]] <- factor(data[[factor_var]], levels = factor_levels)
  } else {
    data[[factor_var]] <- as.factor(data[[factor_var]])
  }

  # 2. Statistical Models
  formula_str <- as.formula(paste(numeric_var, "~", factor_var))
  aov_model <- aov(formula_str, data = data)
  aov_summary <- summary(aov_model)
  tukey_res <- stats::TukeyHSD(aov_model)

  # 3. Estimated Marginal Means & Letters
  emmean <- emmeans::emmeans(aov_model, specs = as.formula(paste("~", factor_var)))
  emmean_cld <- multcomp::cld(emmean, Letters = letters)
  emmean_cld$.group <- trimws(emmean_cld$.group)

  # 4. Summary Statistics for Error Bars & Letter Placement
  sum_data <- data %>%
    dplyr::group_by(.data[[factor_var]]) %>%
    dplyr::summarise(
      mean_val = mean(.data[[numeric_var]], na.rm = TRUE),
      sd_val = sd(.data[[numeric_var]], na.rm = TRUE),
      se_val = sd_val / sqrt(dplyr::n())
    ) %>%
    dplyr::mutate(
      error_margin = if (error_type == "sd") sd_val else se_val,
      upper_limit = mean_val + error_margin
    )

  emmean_cld <- merge(emmean_cld, sum_data, by = factor_var)

  # 5. Base Plot Initialization
  p <- ggplot(data, aes(x = .data[[factor_var]], y = .data[[numeric_var]], fill = .data[[factor_var]])) +
    theme_bw() +
    labs(x = factor_var, y = numeric_var, title = "One-Way ANOVA") +
    scale_fill_brewer(palette = color_palette) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  # 6. Plot Type Routing
  if (plot_type == "boxplot") {
    p <- p + geom_boxplot(outlier.shape = NA, alpha = 0.7)
    if (add_jitter) p <- p + geom_jitter(color = "black", size = 1.5, width = 0.2, alpha = 0.5)
    if (show_mean) p <- p + stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = mean_color)

  } else if (plot_type == "barplot") {
    p <- ggplot(sum_data, aes(x = .data[[factor_var]], y = mean_val, fill = .data[[factor_var]])) +
      geom_col(alpha = 0.8, color = "black") +
      geom_errorbar(aes(ymin = mean_val - error_margin, ymax = upper_limit), width = 0.2) +
      theme_bw() +
      labs(x = factor_var, y = numeric_var, title = paste("One-Way ANOVA (Mean \u00B1", toupper(error_type), ")")) +
      scale_fill_brewer(palette = color_palette) +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  } else if (plot_type == "pointrange") {
    p <- ggplot(sum_data, aes(x = .data[[factor_var]], y = mean_val, color = .data[[factor_var]])) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_val - error_margin, ymax = upper_limit), width = 0.2) +
      theme_bw() +
      labs(x = factor_var, y = numeric_var, title = paste("One-Way ANOVA (Mean \u00B1", toupper(error_type), ")")) +
      scale_color_brewer(palette = color_palette) +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  }

  # 7. Significance Display
  if (sig_display == "letters") {
    y_nudge <- max(data[[numeric_var]], na.rm = TRUE) * 0.05
    if (plot_type == "boxplot") {
      max_y <- max(data[[numeric_var]], na.rm = TRUE)
      p <- p + geom_text(data = emmean_cld, aes(x = .data[[factor_var]], y = max_y + y_nudge, label = .group), inherit.aes = FALSE, size = 5)
    } else {
      p <- p + geom_text(data = emmean_cld, aes(x = .data[[factor_var]], y = upper_limit + y_nudge, label = .group), inherit.aes = FALSE, size = 5)
    }
  } else if (sig_display == "stars") {
    p <- p + ggpubr::geom_pwc(data = data, aes(x = .data[[factor_var]], y = .data[[numeric_var]]),
                              method = "tukey_hsd", label = "p.adj.signif", hide.ns = TRUE, inherit.aes = FALSE,
                              step.increase = 0.1, vjust = 0.5)
  }

  # 8. Apply Custom Y-Axis Limits (if provided)
  if (!is.null(y_limits)) {
    if (length(y_limits) == 2 && is.numeric(y_limits)) {
      p <- p + coord_cartesian(ylim = y_limits)
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

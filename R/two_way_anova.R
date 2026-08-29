#' Two-Way Analysis of Variance (ANOVA) with Interaction Visualization
#'
#' This function performs a two-way ANOVA with interaction effects, calculates
#' Tukey's HSD post-hoc tests, extracts Estimated Marginal Means (EMMs), and
#' generates publication-ready interaction plots with compact letter displays.
#'
#' @param data A data frame in long format.
#' @param factor1_var Character; the name of the first categorical independent variable (x-axis).
#' @param factor2_var Character; the name of the second categorical independent variable (groups/fill).
#' @param numeric_var Character; the name of the continuous response variable.
#' @param factor1_levels Optional character vector; custom order for the first factor levels.
#' @param factor2_levels Optional character vector; custom order for the second factor levels.
#' @param plot_type Character; type of interaction plot: "boxplot", "barplot", or "pointrange". Default is "boxplot".
#' @param error_type Character; error bar type for bar/pointrange plots: "se" or "sd". Default is "se".
#' @param y_limits Optional numeric vector of length 2; explicitly sets the Y-axis limits.
#' @param add_jitter Logical; if TRUE, adds jittered points to the boxplot. Default is TRUE.
#' @param show_mean Logical; if TRUE, displays a mean point inside the boxplot. Default is TRUE.
#' @param mean_color Character; color for the mean point in the boxplot. Default is "darkred".
#' @param color_palette Character; name of a color palette from RColorBrewer. Default is "Set1".
#'
#' @return A list containing the ANOVA summary, Tukey HSD results, EMMs, summary statistics table, and the ggplot object.
#'
#' @import ggplot2 emmeans multcomp RColorBrewer dplyr
#' @importFrom stats aov as.formula sd
#' @importFrom dplyr %>%
#' @export
two_way_anova <- function(data,
                          factor1_var,
                          factor2_var,
                          numeric_var,
                          factor1_levels = NULL,
                          factor2_levels = NULL,
                          plot_type = "boxplot",
                          error_type = "se",
                          y_limits = NULL,
                          add_jitter = TRUE,
                          show_mean = TRUE,
                          mean_color = "darkred",
                          color_palette = "Set1") {

  # 1. Input Validation & Factor Ordering
  required_cols <- c(factor1_var, factor2_var, numeric_var)
  if (!all(required_cols %in% names(data))) stop("Specified columns not found in data.")
  if (!is.numeric(data[[numeric_var]])) stop(paste(numeric_var, "must be numeric."))

  if (!is.null(factor1_levels)) {
    data[[factor1_var]] <- factor(data[[factor1_var]], levels = factor1_levels)
  } else {
    data[[factor1_var]] <- as.factor(data[[factor1_var]])
  }

  if (!is.null(factor2_levels)) {
    data[[factor2_var]] <- factor(data[[factor2_var]], levels = factor2_levels)
  } else {
    data[[factor2_var]] <- as.factor(data[[factor2_var]])
  }

  # 2. Statistical Models
  formula_str <- as.formula(paste(numeric_var, "~", factor1_var, "*", factor2_var))
  aov_model <- aov(formula_str, data = data)
  aov_summary <- summary(aov_model)
  tukey_res <- TukeyHSD(aov_model)

  # 3. Estimated Marginal Means & Letters for Interaction
  emmean <- emmeans::emmeans(aov_model, specs = as.formula(paste("~", factor1_var, "*", factor2_var)))
  emmean_cld <- multcomp::cld(emmean, Letters = letters)
  emmean_cld$.group <- trimws(emmean_cld$.group)

  # 4. Summary Statistics for Error Bars
  sum_data <- data %>%
    dplyr::group_by(.data[[factor1_var]], .data[[factor2_var]]) %>%
    dplyr::summarise(
      mean_val = mean(.data[[numeric_var]], na.rm = TRUE),
      sd_val = sd(.data[[numeric_var]], na.rm = TRUE),
      se_val = sd_val / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      error_margin = if (error_type == "sd") sd_val else se_val,
      upper_limit = mean_val + error_margin
    )

  emmean_cld <- merge(emmean_cld, sum_data, by = c(factor1_var, factor2_var))

  # 5. Base Plot Initialization
  p <- ggplot(data, aes(x = .data[[factor1_var]], y = .data[[numeric_var]], fill = .data[[factor2_var]])) +
    theme_bw() +
    labs(x = factor1_var, y = numeric_var, fill = factor2_var, title = "Two-Way ANOVA Interaction") +
    scale_fill_brewer(palette = color_palette) +
    theme(plot.title = element_text(hjust = 0.5))

  # 6. Plot Type Routing
  if (plot_type == "boxplot") {
    p <- p + geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8))
    if (add_jitter) {
      p <- p + geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
                          color = "black", size = 1.2, alpha = 0.5, show.legend = FALSE)
    }
    if (show_mean) {
      p <- p + stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
                            color = mean_color, position = position_dodge(0.8), show.legend = FALSE)
    }

  } else if (plot_type == "barplot") {
    p <- ggplot(sum_data, aes(x = .data[[factor1_var]], y = mean_val, fill = .data[[factor2_var]])) +
      geom_col(alpha = 0.8, color = "black", position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = mean_val - error_margin, ymax = upper_limit),
                    width = 0.2, position = position_dodge(0.8)) +
      theme_bw() +
      labs(x = factor1_var, y = numeric_var, fill = factor2_var, title = paste("Two-Way ANOVA (Mean \u00B1", toupper(error_type), ")")) +
      scale_fill_brewer(palette = color_palette) +
      theme(plot.title = element_text(hjust = 0.5))

  } else if (plot_type == "pointrange") {
    p <- ggplot(sum_data, aes(x = .data[[factor1_var]], y = mean_val, color = .data[[factor2_var]])) +
      geom_point(size = 3, position = position_dodge(0.4)) +
      geom_errorbar(aes(ymin = mean_val - error_margin, ymax = upper_limit),
                    width = 0.2, position = position_dodge(0.4)) +
      theme_bw() +
      labs(x = factor1_var, y = numeric_var, color = factor2_var, title = paste("Two-Way ANOVA (Mean \u00B1", toupper(error_type), ")")) +
      scale_color_brewer(palette = color_palette) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  # 7. Significance Display (Compact Letter Display)
  y_nudge <- max(data[[numeric_var]], na.rm = TRUE) * 0.05
  if (plot_type == "boxplot") {
    max_y <- max(data[[numeric_var]], na.rm = TRUE)
    p <- p + geom_text(data = emmean_cld, aes(x = .data[[factor1_var]], y = max_y + y_nudge, label = .group, group = .data[[factor2_var]]),
                       position = position_dodge(0.8), inherit.aes = FALSE, size = 4)
  } else {
    p <- p + geom_text(data = emmean_cld, aes(x = .data[[factor1_var]], y = upper_limit + y_nudge, label = .group, group = .data[[factor2_var]]),
                       position = position_dodge(0.8), inherit.aes = FALSE, size = 4)
  }

  # 8. Apply Custom Y-Axis Limits
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

